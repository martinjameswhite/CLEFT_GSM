#include	<cmath>
#include	<complex>
#include	<cstdlib>
#include	<vector>
#include	<fstream>
#include	<sstream>
#include	<iostream>
#include	<iomanip>
#include	<string>
#include	<exception>
#include	"omp.h"
#include	"fftw3.h"




constexpr int Nfft=2048;


// These can be modified later if we convert to MPI.
void	myexit(const int flag) {
  exit(flag);
}
void	myexception(const std::exception& e) {
  std::cout<<"Exception: "<<e.what()<<std::endl;
  myexit(1);
}






class	GaussLegendre {
// Class to hold the abscissae and weights for G-L integration.
// Form:  \int_{-1}^{+1} dx f(x) = \sum w_j f(x_j)
public:
  std::vector<double>	x,w;
  int			N;
  void set(const int NN) {
    // Set the weights and abscissae, which are public variables.
    N = NN;
    try {
      x.resize(N);
      w.resize(N);
    } catch(std::exception& e) {myexception(e);}
    // First find the abscissae...this involves finding the roots
    // of the Legendre polynomials.  We use Newton-Raphson starting
    // from an analytic guess.
    int N2 = (N+1)/2;	// Weights/roots are symmetric, do half ...
    for (int i=0; i<N2; ++i) {
      // Find root by starting "close" and using N-R.
      double zz,dpdz,z=cos(M_PI*(i+0.75)/(N+0.5));
      do {
        double p1=1,p2=0;
        for (int j=0; j<N; ++j) {
          double p3 = p2;
          p2 = p1;
          p1 = ((2*j+1.)*z*p2-j*p3)/(j+1.);
        }
        // Now p1 is P_n and p2 is P_{n-1}.  Compute dP/dz
        dpdz = N*(z*p1-p2)/(z*z-1);
        // and use N-R to update:
        zz = z;
        z  = zz-p1/dpdz;
      } while(fabs(z-zz)>1e-15);
      x[  i  ] = -z;   w[  i  ] = 2.0/(1-z*z)/dpdz/dpdz;
      x[N-i-1] =  z;   w[N-i-1] = 2.0/(1-z*z)/dpdz/dpdz;
    }
  }
  GaussLegendre() {}
  GaussLegendre(const int N) {set(N);}
};






class	Spline {
private:
  std::vector<double> xa,ya,y2;
  int	nn;
public:
  Spline(const std::vector<double>& x, const std::vector<double>& y) {
    double p,qn,sig,un;
    if (x.size() != y.size()) {
      std::cout << "x and y must have the same dimensions." << std::endl;
      myexit(1);
    }
    nn = x.size();
    if (x[0]>=x[nn-1]) {
      std::cout << "x must be monotonically increasing in spLine." << std::endl;
      myexit(1);
    }
    xa = x;
    ya = y;
    try {y2.resize(nn);} catch(std::exception& e) {myexception(e);}
    std::vector<double> u(nn);
    y2[0]=u[0]=0.0;
    for (int i=1;i<nn-1;++i) {
      sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p=sig*y2[i-1]+2.0;
      y2[i]=(sig-1.0)/p;
      u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
      u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    qn=un=0.0;
    y2[nn-1]=(un-qn*u[nn-2])/(qn*y2[nn-2]+1.0);
    for (int k=nn-2;k>=0;--k)
      y2[k]=y2[k]*y2[k+1]+u[k];
  }
  double operator() (const double x) {
    int	k,klo,khi;
    if (x<xa[0] || x>xa[nn-1]) {
      std::cout << "x out of range in Spline." << std::endl;
      myexit(1);
    }
    klo=0;
    khi=nn-1;
    while (khi-klo > 1) {	// Bisection search.
      k=(khi+klo) >> 1;
      if (xa[k] > x) khi=k; else klo=k;
    }
    double h=xa[khi]-xa[klo];
    if (h<=0.0) {std::cout << "h==0 in bisection." << std::endl; myexit(1); }
    double a=(xa[khi]-x)/h;
    double b=(x-xa[klo])/h;
    return(a*ya[klo]+b*ya[khi]+
      ((a*a*a-a)*y2[klo]+(b*b*b-b)*y2[khi])*(h*h)/6.0);
  }
  double xmin() {
    return(xa[0]);
  }
  double xmax() {
    return(xa[nn-1]);
  }
};







class LanczosGamma {
// Only one real method, to return the Lanzos approximation to ln[Gamma]
// for complex arguments.  The coefficients are taken from the n=15
// table of https://mrob.com/pub/ries/lanczos-gamma.html .
private:
    const double         LG_g = 4.7421875;
    const double ln_sqrt_2_pi = 0.91893853320467274178;
    const double         g_pi = 3.14159265358979323846;
    const std::vector<double> lct={0.99999999999999709182,
         57.156235665862923517,
        -59.597960355475491248,
         14.136097974741747174,
        -0.49191381609762019978,
        0.33994649984811888699e-4,
        0.46523628927048575665e-4,
        -0.98374475304879564677e-4,
        0.15808870322491248884e-3,
        -0.21026444172410488319e-3,
        0.21743961811521264320e-3,
        -0.16431810653676389022e-3,
        0.84418223983852743293e-4,
        -0.26190838401581408670e-4,
        0.36899182659531622704e-5};
public:
  std::complex<double> ln_gamma(const std::complex<double> z) const {
    if (z.real() < 0.5) {
      // Use Euler's reflection formula:
      // Gamma(z) = Pi / [Sin[Pi*z] * Gamma[1-z]];
      return(log(g_pi / sin(g_pi * z)) - ln_gamma(1.0 - z));
    }
    std::complex<double> zz   = z - 1.0;
    std::complex<double> base = zz + LG_g + 0.5;
    std::complex<double> sum  = {0,0};
    for(int i=lct.size()-1; i>=1; --i)
      sum += lct[i] / (zz + ((double) i));
    sum += lct[0];
    // Gamma[z] = Sqrt(2*Pi) * sum * base^[z + 0.5] / E^base
    std::complex<double> val=ln_sqrt_2_pi+log(sum)-base+log(base)*(zz+0.5);
    return(val);
  }
};





template <int Nfft> class FFTlog {
private:
  const LanczosGamma lg;
  const int L;
  const double sqrtpi=1.7724538509055160272981674833411451827975494561223865;
  std::vector<double> fk;
  std::vector< std::vector<double> >   ys,us;
  double	Delta;
  fftw_plan	fplan,iplan;
  //
  std::complex<double> UK(const std::complex<double>& nu,
                          const std::complex<double>& z) const {
    std::complex<double> ret;
    ret = log(2.0)*(z-2.0)+lg.ln_gamma(0.5*(nu+z))-lg.ln_gamma(0.5*(3.0+nu-z));
    return(sqrtpi*exp(ret));
  }
public:
  // These are public so we can write directly into fq to save a little time.
  std::vector<double> qs,fq,ks;
  FFTlog(const std::vector<double>& qq, const int Lmax=15):L(Lmax),lg() {
    // Initialize with the values of the untransformed coordinate, q,
    // and an Lmax value.
    if (qq.size()!=Nfft) {
      std::cout<<"FFTlog vector q.size()="<<qq.size()
	       <<" != Nfft="<<Nfft<<std::endl;
      myexit(1);
    }
    qs    = qq;
    Delta = log(qs[Nfft-1]/qs[0])/(Nfft-1);
    // Allocate memory for FFT arrays.
    ks.resize(Nfft);
    fq.resize(Nfft);
    fk.resize(Nfft+2);
    fftw_complex *fout=(fftw_complex *)&fk[0];
    // Set up padding.
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
    fplan = fftw_plan_dft_r2c_1d(Nfft,&fq[0],fout,FFTW_ESTIMATE);
    iplan = fftw_plan_dft_c2r_1d(Nfft,fout,&fq[0],FFTW_ESTIMATE);
  }
  ~FFTlog() {
    fftw_destroy_plan(fplan);
    fftw_destroy_plan(iplan);
    fftw_cleanup_threads();
  }
  void init() {
    // Set up the L arrays holding the Y and Mellin transforms u_m.
    // We store the u_m, which are complex, as real-imag-real-imag-real...
    // since this is how the FFT will access them.
    ys.clear(); us.clear();
    for (int ell=0; ell<L; ++ell) {
      std::vector<double> yy(Nfft),um(Nfft+2);
      double qval = (1.5<ell)?0:1.5-ell;
      std::complex<double> qrot(qval,M_PI/Delta);
      std::complex<double> uofk=UK(ell+0.0,qrot);
      double lnxy = Delta/M_PI*atan2(uofk.imag(),uofk.real());
      for (int i=0; i<Nfft; ++i)
        yy[i]=exp(lnxy-Delta)*qs[i]/(qs[0]*qs[Nfft-1]);
      ys.push_back(yy);
      const double fact=2*M_PI/Nfft/Delta;
      for (int j=0; j<Nfft/2+1; ++j) {
        std::complex<double> tmp1(qval,j*fact);
        std::complex<double> tmp2(0.0,-j*fact*lnxy);
        std::complex<double> tmp3=UK(ell+0.0,tmp1)*exp(tmp2);
        um[2*j+0]=tmp3.real();
        um[2*j+1]=tmp3.imag();
      }
      um[Nfft+1] = 0;  // Enforce reality condition.
      us.push_back(um);
    }
  }
  void sph(const int nu, const bool k2q) {
    // The main method: does the spherical Hankel transform.  The
    // k2q flag divides by (2\pi^2) if going from k->q.
    const double qval=(1.5<nu)?0:1.5-nu;
    std::vector<double>& um = us[nu];
    // Multiply the f(q) by qs^{3-qval}
    for (int i=0; i<Nfft; ++i)
      fq[i] *= pow(qs[i],3.0-qval);
    // Do the forward transform, from fq->fk.
    fftw_execute(fplan);
    // Do the convolution:
    for (int i=0; i<Nfft/2+1; ++i) {
      double tmpr,tmpi;
      tmpr = (fk[2*i+0]*um[2*i+0]-fk[2*i+1]*um[2*i+1]);
      tmpi =-(fk[2*i+0]*um[2*i+1]+fk[2*i+1]*um[2*i+0]);
      fk[2*i+0] = tmpr;
      fk[2*i+1] = tmpi;
    }
    // Transform back.
    fftw_execute(iplan);
    // Rescale the output.
    const double norm = (k2q)?1.0/Nfft/2/M_PI/M_PI:1.0/Nfft;
    std::vector<double>& yy = ys[nu];
    for (int i=0; i<Nfft; ++i) {
      ks[i]  = yy[i];
      fq[i] *= pow(yy[i],-qval) * norm;
    }
  }
};




class	Zeldovich {
// Computes the power spectrum within the Zel'dovich approximation.
protected:
  const int NkTable=8192;
  std::vector<double>	kLin,pLin;
  std::vector<double>	kval,Pval;
  double dkinv;
  bool   doneQfuncs;
  void readPowerSpectrum(const char fname[]) {
    // Load the linear power spectrum from a file, and expand it out
    // if necessary.  This is stored log-log.
    kLin.clear();  pLin.clear();
    std::ifstream fs(fname);
    if (!fs) {
      std::cout<<"Unable to open "<<fname<<" for reading."<<std::endl;
      myexit(1);
    }
    std::string buf;
    do {	// Skip any preamble comment lines.
      getline(fs,buf);
    } while(!fs.eof() && buf[0]=='#');
    while (!fs.eof()) {
      double kval,pval;
      std::istringstream ss(buf);
      ss >> kval >> pval;
      try {
        kLin.push_back(log(kval));
        pLin.push_back(log(pval));
      } catch(std::exception& e) {myexception(e);}
      getline(fs,buf);
    }
    fs.close();
    // Now resample to an "even" log-k spacing.
    Spline ss(kLin,pLin);
    kLin.clear(); pLin.clear();
    for (int i=0; i<NkTable; ++i) {
      double x = ss.xmin() + (i+0.5)*(ss.xmax()-ss.xmin())/NkTable;
      try {
        kLin.push_back(x);
        pLin.push_back(ss(x));
      } catch(std::exception& e) {myexception(e);}
    }
    // and set up dkinv so we can interpolate quickly.
    dkinv = NkTable/(kLin[NkTable-1]-kLin[0]);
  }
  void makePkarray() {
    const int    Npad=256;
    const double lnKmin=log(1e-5);
    const double lnKmax=log(1e5 );
    const double dlnk=(lnKmax-lnKmin)/(Nfft-2*Npad);
#pragma omp parallel for shared(kval)
    for (int i=0; i<Nfft; ++i)
      kval[i] = exp( lnKmin+(i-Npad)*dlnk );
    // Apodize the P(k) to lie within the window:
#pragma omp parallel for shared(kval,Pval)
    for (int i=0; i<Nfft; ++i)
      Pval[i] = linearPk(kval[i]) *
                0.5*(1.0+tanh(0.125*(i-Npad))) *
                0.5*(1.0+tanh(0.125*(Nfft-Npad-i)));
  }
public:
  std::vector<double>	qqLin,xiLin,UULin,XXLin,YYLin,YqLin;
  std::vector<double>	XsLin,YsLin,VVLin,zetaLin,chiLin;
  double sigma2,shear_sig2;
  Zeldovich(const char fname[]) {
    doneQfuncs = false;
    kval.resize(Nfft);
    Pval.resize(Nfft);
    readPowerSpectrum(fname);
    makePkarray();
  }
  double linearPk(const double kk) const {
    // Returns the interpolated, linear theory power spectrum or a
    // power-law extrapolation (if outside the bounds).
    double xx=log(kk);
    int    jj=(int)((xx-kLin[0])*dkinv);
    if (jj<=0) jj=0;
    if (jj>=pLin.size()-2) jj=pLin.size()-2;
    double pk = exp(pLin[jj]+(xx-kLin[jj])*
                   (pLin[jj+1]-pLin[jj])/(kLin[jj+1]-kLin[jj]));
    return(pk);
  }
  void calcQfuncs() {
    // Computes the q-dependent functions, X, Y, U and xi.
    // Set up a regular q grid:
    qqLin.resize(Nfft);
    xiLin.resize(Nfft);
    XXLin.resize(Nfft);
    YYLin.resize(Nfft);
    YqLin.resize(Nfft);
    const double lnqmin=log(1e-3);
    const double lnqmax=log(3e4);
#pragma omp parallel for shared(qqLin)
    for (int i=0; i<Nfft; ++i)
      qqLin[i] = exp( lnqmin+(i+0.5)*(lnqmax-lnqmin)/Nfft );
    // Set up the FFTlog class.
    const int Lmax=5;
    FFTlog<Nfft> fftlog(kval,Lmax);
    fftlog.init();
    // Set up some useful powers of kk.
    std::vector<double> kp2(Nfft),km2(Nfft);
#pragma omp parallel for shared(kval,kp2,km2)
    for (int i=0; i<Nfft; ++i) {
      kp2[i] =     kval[i]*kval[i];
      km2[i] = 1.0/kval[i]/kval[i];
    }
    // First compute the xi_ell^n:
    std::vector<double>	xi0m0(Nfft),xi0m2(Nfft),xi1m1(Nfft),xi2m2(Nfft);
    // extra xi_ell^n for shear terms:
    std::vector<double>	xi2p0(Nfft),xi4p0(Nfft),xi1p1(Nfft);
    std::vector<double>	xi3p1(Nfft),xi3m1(Nfft),xi0p2(Nfft);
    // We could in principle do this with a std::map to keep all of
    // the xi_ell_n arrays neatly organized, but this seemed easier
    // at the time:
    // ell=0; n=0:
    if (true) {
      fftlog.fq = Pval;
      fftlog.sph(0,true);
      // and sample onto a regular grid.
      Spline ss(fftlog.ks,fftlog.fq);
#pragma omp parallel for shared(qqLin,ss,xi0m0)
      for (int i=0; i<Nfft; ++i)
        xi0m0[i] = ss(qqLin[i]);
    }
    // ell=0; n=-2:
    if (true) {
      for (int i=0; i<Nfft; ++i)
        fftlog.fq[i] = Pval[i] * km2[i];
      fftlog.sph(0,true);
      Spline ss(fftlog.ks,fftlog.fq);
#pragma omp parallel for shared(qqLin,ss,xi0m2)
      for (int i=0; i<Nfft; ++i)
        xi0m2[i] = ss(qqLin[i]);
    }
    // ell=1; n=-1:
    if (true) {
      for (int i=0; i<Nfft; ++i)
        fftlog.fq[i] = Pval[i] / kval[i];
      fftlog.sph(1,true);
      Spline ss(fftlog.ks,fftlog.fq);
#pragma omp parallel for shared(qqLin,ss,xi1m1)
      for (int i=0; i<Nfft; ++i)
        xi1m1[i] = ss(qqLin[i]);
    }
    // ell=2; n=-2:
    if (true) {
      for (int i=0; i<Nfft; ++i)
        fftlog.fq[i] = Pval[i] * km2[i];
      fftlog.sph(2,true);
      Spline ss(fftlog.ks,fftlog.fq);
#pragma omp parallel for shared(qqLin,ss,xi2m2)
      for (int i=0; i<Nfft; ++i)
        xi2m2[i] = ss(qqLin[i]);
    }
    // Now the shear terms...
    // ell=2; n=0:
    if (true) {
      fftlog.fq = Pval;
      fftlog.sph(2,true);
      Spline ss(fftlog.ks,fftlog.fq);
#pragma omp parallel for shared(qqLin,ss,xi2p0)
      for (int i=0; i<Nfft; ++i)
        xi2p0[i] = ss(qqLin[i]);
    }
    // ell=4; n=0:
    if (true) {
      fftlog.fq = Pval;
      fftlog.sph(4,true);
      Spline ss(fftlog.ks,fftlog.fq);
#pragma omp parallel for shared(qqLin,ss,xi4p0)
      for (int i=0; i<Nfft; ++i)
        xi4p0[i] = ss(qqLin[i]);
    }
    // ell=1; n=+1:
    if (true) {
      for (int i=0; i<Nfft; ++i)
        fftlog.fq[i] = Pval[i] * kval[i];
      fftlog.sph(1,true);
      Spline ss(fftlog.ks,fftlog.fq);
#pragma omp parallel for shared(qqLin,ss,xi1p1)
      for (int i=0; i<Nfft; ++i)
        xi1p1[i] = ss(qqLin[i]);
    }
    // ell=3; n=+1:
    if (true) {
      for (int i=0; i<Nfft; ++i)
        fftlog.fq[i] = Pval[i] * kval[i];
      fftlog.sph(3,true);
      Spline ss(fftlog.ks,fftlog.fq);
#pragma omp parallel for shared(qqLin,ss,xi3p1)
      for (int i=0; i<Nfft; ++i)
        xi3p1[i] = ss(qqLin[i]);
    }
    // ell=3; n=-1:
    if (true) {
      for (int i=0; i<Nfft; ++i)
        fftlog.fq[i] = Pval[i] / kval[i];
      fftlog.sph(3,true);
      Spline ss(fftlog.ks,fftlog.fq);
#pragma omp parallel for shared(qqLin,ss,xi3m1)
      for (int i=0; i<Nfft; ++i)
        xi3m1[i] = ss(qqLin[i]);
    }
    // ell=0; n=+2:
    if (true) {
      for (int i=0; i<Nfft; ++i)
        fftlog.fq[i] = Pval[i] * kp2[i];
      fftlog.sph(0,true);
      Spline ss(fftlog.ks,fftlog.fq);
#pragma omp parallel for shared(qqLin,ss,xi0p2)
      for (int i=0; i<Nfft; ++i)
        xi0p2[i] = ss(qqLin[i]);
    }
    // Now set up the "q" functions.
    XXLin.resize(Nfft);
    YYLin.resize(Nfft);
    xiLin.resize(Nfft);
    UULin.resize(Nfft);
    YqLin.resize(Nfft);
#pragma omp parallel for shared(XXLin,YYLin,xiLin,UULin,YqLin,xi0m2,xi2m2,xi0m0,xi1m1)
    for (int i=0; i<Nfft; ++i) {
      XXLin[i] = 2./3.*( xi0m2[0]-xi0m2[i]-xi2m2[i]);
      YYLin[i] = 2*xi2m2[i];
      xiLin[i] =   xi0m0[i];
      UULin[i] =  -xi1m1[i];
      YqLin[i] = YYLin[i]/qqLin[i];
    }
    sigma2 = XXLin[Nfft-1]+YYLin[Nfft-1];
    // Do the shear functions.
      XsLin.resize(Nfft);
      YsLin.resize(Nfft);
      VVLin.resize(Nfft);
     chiLin.resize(Nfft);
    zetaLin.resize(Nfft);
    for (int i=0; i<Nfft; ++i) {
      double  J2 = 2./15.*xi1m1[i]-1./5.*xi3m1[i];
      double  J3 =-1./ 5.*xi1m1[i]-1./5.*xi3m1[i];
      double  J4 = xi3m1[i];
        XsLin[i] = 4*J3*J3;
        YsLin[i] = 6*J2*J2+8*J2*J3+4*J2*J4+4*J3*J3+8*J3*J4+2*J4*J4;
        VVLin[i] = 4.*xi2p0[i]*J2;
       chiLin[i] = 4./3.*xi2p0[i]*xi2p0[i];
      zetaLin[i] = 2*(4./45.*xi0m0[i]*xi0m0[i]+8./63.*xi2p0[i]*xi2p0[i]+
		      8./35.*xi4p0[i]*xi4p0[i]);
    }
    shear_sig2 = XsLin[Nfft-1]+YsLin[Nfft-1];
    doneQfuncs = true;
  }
  void printQfuncs() {
    if (!doneQfuncs) calcQfuncs();
    std::cout<<"# Linear theory, q-dependent functions."<<std::endl;
    std::cout<<"# Sigma2="<<sigma2<<std::endl;
    std::cout<<"# "
	     <<std::setw(10)<<"q"
	     <<std::setw(12)<<"xi"
	     <<std::setw(12)<<"X"
	     <<std::setw(12)<<"Y"
	     <<std::setw(12)<<"U"
	     <<std::endl;
    for (int i=0; i<Nfft; ++i)
        std::cout<<std::scientific<<std::setw(12)<<std::setprecision(4)
                 <<qqLin[i]
                 <<std::scientific<<std::setw(12)<<std::setprecision(4)
                 <<xiLin[i]
                 <<std::scientific<<std::setw(12)<<std::setprecision(4)
                 <<XXLin[i]
                 <<std::scientific<<std::setw(12)<<std::setprecision(4)
                 <<YYLin[i]
                 <<std::scientific<<std::setw(12)<<std::setprecision(4)
                 <<UULin[i]
                 <<std::endl;
  }
  double power(const double kk,
               const double b1=0, const double b2=0,
	       const double bs=0, const double alpha=0) const {
    // Returns the Zeldovich power spectrum at kk.
    if (!doneQfuncs) return(0.0); // Since method is const, can't make Qfuncs.
    const double k2=kk*kk;
    const int Lmax=8;
    FFTlog<Nfft> fftlog(qqLin,Lmax);
    fftlog.init();
    // Set up an apodization function.
    std::vector<double> ap(Nfft);
    const int Npad=256;
#pragma omp parallel for shared(ap)
    for (int i=0; i<Nfft; ++i)
      ap[i] = 0.5*(1.0+tanh(0.125*(i-Npad))) *
              0.5*(1.0+tanh(0.125*(Nfft-Npad-i)));
    // Now compute P(k) as a sum of Hankel transforms:
    double ret=0,pk=0;
    for (int ell=0; ell<Lmax; ++ell) {
      if (ell==0) {
        for (int i=0; i<Nfft; ++i) {
          double U2  = UULin[i]*UULin[i];
          double tmp = -0.5*k2*(XXLin[i]+YYLin[i]-sigma2);
          fftlog.fq[i] =(expm1(tmp) + exp(tmp)*(
                         b1* 1*( 0 ) +
                         b2* 1*( -k2*U2 ) +
                         b1*b1*( xiLin[i]-k2*U2 ) +
                         b1*b2*( 0 ) +
                         b2*b2*( 0.5*xiLin[i]*xiLin[i] ) +
			 bs* 1*( -k2*(XsLin[i]+YsLin[i]) ) +
			 b1*bs*( 0 ) +
			 b2*bs*(  chiLin[i] ) +
			 bs*bs*( zetaLin[i] ) 
			 ) ) * ap[i] + bs*1*k2*shear_sig2*ap[i];
        }
      }
      else {
        for (int i=0; i<Nfft; ++i) {
          double U2  = UULin[i]*UULin[i];
          double tmp = -0.5*k2*(XXLin[i]+YYLin[i]-sigma2);
	  double fac = 1-2*ell/k2/YYLin[i];
          fftlog.fq[i] = exp(tmp) * pow(YqLin[i],double(ell))*(1 +
                         b1* 1*( -2*UULin[i]/YqLin[i] ) +
                         b2* 1*( (2*ell/YYLin[i]-k2)*U2 ) +
                         b1*b1*( xiLin[i]+(2*ell/YYLin[i]-k2)*U2 ) +
                         b1*b2*( -2*xiLin[i]*UULin[i]/YqLin[i] ) +
                         b2*b2*( 0.5*xiLin[i]*xiLin[i] ) +
			 bs* 1*( -k2*(XsLin[i]+fac*YsLin[i]) ) +
			 b1*bs*( -2*VVLin[i]/YqLin[i] ) +
			 b2*bs*(  chiLin[i] ) +
			 bs*bs*( zetaLin[i] )
			 ) * ap[i];
        }
      }
      fftlog.sph(ell,false);
      // Now linearly interpolate the desired value
      int ilo=0,ihi=Nfft-1;
      while (ihi>ilo+1) {
        int imid=(ilo+ihi)/2;
        if (fftlog.ks[imid]<kk) {
          ilo=imid;
        }
        else {
          ihi=imid;
        }
      }
      pk   = fftlog.fq[ilo]+(kk-fftlog.ks[ilo])*
            (fftlog.fq[ihi]-fftlog.fq[ilo])/(fftlog.ks[ihi]-fftlog.ks[ilo]);
      pk  *= pow(kk,double(ell));
      ret += pk;
    }
    ret *= 4*M_PI;
    // Add in an approximation to the EFT term.
    ret += alpha*(k2*sigma2)*linearPk(kk);
    ret *= exp(-0.5*k2*sigma2);
    return(ret);
  }
};






int	main(int argc, char **argv)
{
  if (argc!=6) {
    std::cout<<"Usage: zeldovich <pkfile> <b1> <b2> <bs> <alpha>"<<std::endl;
    exit(1);
  }
  Zeldovich zel(argv[1]);
  double b1 = atof(argv[2]);
  double b2 = atof(argv[3]);
  double bs = atof(argv[4]);
  double alpha = atof(argv[5]);

  zel.calcQfuncs();
  //zel.printQfuncs();

  // Just do an example here.
  std::cout<<"# Zeldovich (real-space) power spectrum."<<std::endl;
  std::cout<<"# b1="<<b1<<", b2="<<b2
           <<", bs="<<bs<<", alpha="<<alpha<<std::endl;
  const int Nk=100;
  for (int i=0; i<Nk; ++i) {
    double kk = exp( log(1e-2)+(i+0.5)*log(3.0/1e-2)/Nk );
    std::cout<<std::scientific<<std::setw(15)<<std::setprecision(5)<<kk
             <<std::scientific<<std::setw(15)<<std::setprecision(5)
             <<zel.power(kk,b1,b2,bs,alpha)
	     <<std::endl;
  }

  return(0);
}
