#include	<cstdlib>
#include	<fstream>
#include	<sstream>
#include	<iostream>
#include	<iomanip>
#include	<cmath>
#include	<vector>
#include	<string>
#include	<exception>
#include	<stdexcept>
#include	<omp.h>


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
    const int N2 = (N+1)/2;	// Weights/roots are symmetric, do half ...
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
  Spline() {nn=0;}
  Spline(const std::vector<double>& x, const std::vector<double>& y) {
    init(x,y);
  }
  void init(const std::vector<double>& x, const std::vector<double>& y) {
    double p,qn,sig,un;
    if (x.size() != y.size()) {
      std::cout << "x and y must have the same dimensions." << std::endl;
      myexit(1);
    }
    nn = x.size();
    if (x[0]>=x[nn-1]) {
      std::cout << "x must be monotonically increasing in Spline." << std::endl;
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
  double operator() (const double x) const {
    int	k,klo,khi;
    if (x<xa[0] || x>xa[nn-1]) {
      std::cout<<"x out of range ["<<xa[0]<<","<<xa[nn-1]<<"] in Spline."
               <<std::endl;
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
  double xmin() const {
    return(xa[0]);
  }
  double xmax() const {
    return(xa[nn-1]);
  }
};




class	LPT {
private:
  std::vector<double>	kLin,pLin;
  double		dkinv;
  double linInterp(const double kk) const {
    // Returns a linear interpolant to P(k).
    double lnk = log(kk);
    if (lnk<=kLin[0] || lnk>=kLin[kLin.size()-1]) {return(0);}
    int jj = (int)((lnk-kLin[0])*dkinv);
    if (jj>=pLin.size()-2) jj=pLin.size()-2;
    double pk  = exp(pLin[jj]+(lnk-kLin[jj])*
                    (pLin[jj+1]-pLin[jj])/(kLin[jj+1]-kLin[jj]));
    return(pk);
  }
  std::vector<double> Qdefs(const double r, const double x) const {
    // The 13 Q kernels -- to keep the numbering simple
    // have an empty Q[0].
    // We have added a 14th Q, beyond Matsubara, for shear.
    std::vector<double> Q(15);
    double  x2,r2,rx,yy,yy2;
    r2    = r*r;
    x2    = x*x;
    rx    = r*x;
    yy    = 1.0+r2-2*rx;
    yy2   = yy*yy;
    Q[ 1] = r2 * (1-x2)*(1-x2)/yy2;
    Q[ 2] = (1-x2)*rx*(1-rx)/yy2;
    Q[ 3] = x2*(1-rx)*(1-rx)/yy2;
    Q[ 4] = (1-x2)/yy2;
    Q[ 5] = rx*(1-x2)/yy;
    Q[ 6] = (1-3*rx)*(1-x2)/yy;
    Q[ 7] = x2 * (1-rx)/yy;
    Q[ 8] = r2 * (1-x2)/yy;
    Q[ 9] = rx*(1-rx)/yy;
    Q[10] = 1-x2;
    Q[11] = x2;
    Q[12] = rx;
    Q[13] = r2;
    Q[14] = r2*(x2-1)*(1-2*r2+4*rx-3*x2)/yy2;
    return(Q);
  }
  std::vector<double> Rdefs(const double r) const {
    // The 2 R kernels.
    std::vector<double> R(2);
    const double r2 = r*r;
    if (r<1e-3) {
      R[0] = r2*(16./15.+r2*(-16./35.+r2*16./315.));
      R[1] = r2*( 4./15.+r2*(-12./35.+r2* 4./ 63.));
    }
    else {
      if (fabs(r-1.0)<1e-6) {
        R[0] = 2./3.;
        R[1] = 0;
      }
      else {
        if (r>100.) {
          double x2   = 1.0/r2;
          R[0] = 16./15. + x2*(-16./35.+x2*16./315.);
          R[1] =- 4./15. + x2*( 12./35.-x2* 4./ 63.);
        }
        else {
          R[0] = -(1+r2)*(3-14*r2+3*r2*r2)/(24*r2) +
                  (1-r2)*(1-r2)*(1-r2)*(1-r2)*
                  log(fabs((1+r)/(1-r)))/(16.*r*r2);
          R[1] =  (1-r2)*(3-2*r2+3*r2*r2)/(24*r2) +
                  (r2-1)*(r2-1)*(r2-1)*
                  (1+r2)*log(fabs((1+r)/(1-r)))/(16.*r*r2);
        }
      }
    }
    return(R);
  }
  std::vector<double> Qnx(const double kval, const double rr) const {
    // Does the inner x-integral in Eq. (A39) of Matsubara.
    const int NQ=15;
    std::vector<double> Qret(NQ);
    const int NN=6000;
    const double dx=2.0/NN;
    for (int j=0; j<NN; ++j) {
      double xx = -1.0 + (j+0.5)*dx;
      double kk = kval*sqrt(1.0+rr*rr-2*rr*xx);
      double pk = linInterp(kk);
      std::vector<double> Q = Qdefs(rr,xx);
      for (int i=1; i<NQ; ++i) Qret[i] += Q[i]*pk;
    }
    for (int i=1; i<NQ; ++i) Qret[i] *= dx;
    return(Qret);
  }
public:
  std::vector<double> Qn(const double kk) const {
    // Returns the 13 Q_n at kk.  Does the outer k-integral of Eq. (A39).
    // Add one additional Q_n for shear.
    const int NQ=15;
    std::vector<double> sum(NQ),QQ(NQ),Qret(NQ);
    // Do the integral over r in 2 pieces: first 0-1:
    const int NN=2000;
    const double dr = (1.0-0.0)/NN;
    for (int i=1; i<NQ; ++i) sum[i]=0;
    for (int j=0; j<NN; ++j) {
      double rr = 0.0 + (j+0.5)*dr;
      double pk = linInterp(kk*rr);
      QQ = Qnx(kk,rr);
      for (int i=1; i<NQ; ++i) sum[i] += QQ[i]*pk;
    }
    for (int i=1; i<NQ; ++i) Qret[i] = sum[i]*dr;
    // Now 1-infinity using 1/r as the integration variable.
    for (int i=1; i<NQ; ++i) sum[i]=0;
    for (int j=0; j<NN; ++j) {
      double rr = 0.0 + (j+0.5)*dr;
      double pk = linInterp(kk/rr);
      QQ = Qnx(kk,1.0/rr);
      for (int i=1; i<NQ; ++i) sum[i] += QQ[i]*pk/rr/rr;
    }
    for (int i=1; i<NQ; ++i) {
      Qret[i] += sum[i]*dr;
      Qret[i] *= kk*kk*kk/4./M_PI/M_PI;
    }
    return(Qret);
  }
  std::vector<double> Rn(const double kk) const {
    // Does the integral of Eq. (A40).
    std::vector<double> R(2),sum(2);
    const int NN=5000;
    // Do the integral over r in 2 pieces: first 0-1:
    const double dr = (1.0-0.0)/NN;
    for (int i=0; i<2; ++i) sum[i]=0;
    for (int j=0; j<NN; ++j) {
      double rr = 0.0 + (j+0.5)*dr;
      double pk = linInterp(kk*rr);
      std::vector<double> RR = Rdefs(rr);
      sum[0] += RR[0]*pk;
      sum[1] += RR[1]*pk;
    }
    R[0] = sum[0]*dr;
    R[1] = sum[1]*dr;
    // Now 1-infinity using 1/r as the integration variable.
    sum[0]=sum[1]=0;
    for (int j=0; j<NN; ++j) {
      double rr = 0.0 + (j+0.5)*dr;
      double pk = linInterp(kk/rr);
      std::vector<double> RR = Rdefs(1.0/rr);
      sum[0] += RR[0]*pk/rr/rr;
      sum[1] += RR[1]*pk/rr/rr;
    }
    for (int i=0; i<2; ++i) {
      R[i] += sum[i]*dr;
      R[i] *= kk*kk*kk/4./M_PI/M_PI * linInterp(kk);
    }
    return(R);
  }
  void init(std::vector<double> inkLin, std::vector<double> inpLin) {
    // Simply take a copy of kLin/pLin, which should be in log-log
    // format and equally spaced, and set dkinv.
    try {
      kLin.resize(inkLin.size());
      pLin.resize(inpLin.size());
    } catch(std::exception& e) {myexception(e);}
    for (int i=0; i<kLin.size(); ++i) {
      kLin[i] = inkLin[i];
      pLin[i] = inpLin[i];
    }
    dkinv = kLin.size()/(kLin[kLin.size()-1]-kLin[0]);
  }
};











class	Zeldovich {
// Computes the correlation function within the Zel'dovich approximation.
// Based upon Carlson et al. (2013; MNRAS, 429, 1674): CLPT.
protected:
  // Table lengths for storing functions -- should be even.
  const static int	NqTable=2048,NkTable=2048;
  GaussLegendre		gl;
  std::vector<double>	kLin,pLin,etaPer,etaPar,uVal;
  std::vector<double>	xiLin,dxiLin,n2Lin,J2Lin,J3Lin,J4Lin;
  std::vector<double>	V12Lin,chiLin,zetLin;
  double		qmin,qmax,delta,dkinv,sigma2,eftNorm;
  std::vector<double> sphBess(const double x) {
    // Returns j0(x) and j1(x)/x.
    std::vector<double> jl(2);
    if (fabs(x)<1e-3) {
      double x2=x*x;
      jl[0] = 1.0 + x2*(-1/6.  +x2/120.);
      jl[1] = 1/3.+ x2*(-1./30.+x2/840.);
    }
    else {
      jl[0] = sin(x)/x;
      jl[1] = (jl[0]-cos(x))/x/x;
    }
    return(jl);
  }
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
  double calcSigma2() {
    // Computes sigma^2_eta, the (square of the) 1D dispersion in the
    // displacement field.  Eq. (29) of CLPT paper.  Does the integral
    // in ln(k), assuming kLin is equally spaced in ln(k) and that there
    // are enough points in the array for such an integral.
    double sum=0;
#pragma omp parallel for reduction(+:sum)
    for (int i=1; i<kLin.size(); ++i) {
      int wt = 2+2*(i%2);
      sum   += exp(kLin[i]+pLin[i])*wt;
    }
    sum *= (kLin[2]-kLin[0])/6;
    sum /= 6*M_PI*M_PI;
    return(sum);
  }
  double calcEFTnorm() {
    // It is useful to have some semi-consistent normalization of the EFT
    // terms, so that Aeft can be O(1).  One way of doing this is to
    // multiply everything by Sigma^2/xi_0(0).  But we don't want to do
    // this with different factors during Recon, so we compute a scaling that
    // depends only on Plin, not filtered Plin.
    // By removing the xi_L(0) in the xiContributions below we can just do
    // the Sigma^2 piece here.
    double sum=0;
#pragma omp parallel for reduction(+:sum)
    for (int i=1; i<kLin.size(); ++i) {
      double ap = cos(M_PI/2.*exp(kLin[i]-kLin[kLin.size()-1]));
      int    wt = 2+2*(i%2);
      sum      += exp(kLin[i]+1*pLin[i])*ap*wt;         // Sigma^2
    }
    sum *= (kLin[2]-kLin[0])/6;
    sum /= 6*M_PI*M_PI;
    return(sum);
  }
  std::vector<double> calcQfuncs(const double q) {
    // Returns etaPerp and etaPar, Eqs. (30-31) of CLPT paper.
    // These, up to a factor of f, are the \Psi of Reid&White, Eqs. (9-10).
    // Also, since it is useful, returns U(q) of Eq. (32) as qf[2]
    // the linear xi as qf[3] and its derivative as qf[4].
    const int Nk=kLin.size();
    const double kmax=exp(kLin[Nk-1]);
    int Nint=(int)(8*kmax*q+512);
    if (Nint>=10000) Nint=10000;
    const double hh=(kLin[Nk-1]-kLin[0])/Nint;
    double sum0=0,sum1=0,sum2=0,sum3=0;
#pragma omp parallel for reduction(+:sum0,sum1,sum2,sum3)
    for (int i=1; i<Nint; ++i) {
      double xx = kLin[0]+i*hh;
      double ap = cos(M_PI/2.*exp(xx-kLin[Nk-1]));
      double kk = exp(xx);
      double k2 = kk*kk;
      int    jj = (int)(i*hh*dkinv);
      if (jj>=pLin.size()-2) jj=pLin.size()-2;
      double pk = exp(pLin[jj]+(xx-kLin[jj])*
                     (pLin[jj+1]-pLin[jj])/(kLin[jj+1]-kLin[jj]));
      std::vector<double> jl=sphBess(kk*q);
      int wt= 2+2*(i%2);
      sum0 += kk*pk*(           jl[1])*wt;	// eta_per, Eq. (30)
      sum1 += kk*pk*(jl[0]   -2*jl[1])*wt;	// eta_par, Eq. (31)
      sum2 +=-kk*pk*(k2*      q*jl[1])*wt;	// U,       Eq. (32)
      sum3 += kk*pk*(k2*        jl[0])*wt*ap;	// xi_lin
    }
    std::vector<double> sum(4);
    sum[0] = sum0 * hh/3.0/(2*M_PI*M_PI);
    sum[1] = sum1 * hh/3.0/(2*M_PI*M_PI);
    sum[2] = sum2 * hh/3.0/(2*M_PI*M_PI);
    sum[3] = sum3 * hh/3.0/(2*M_PI*M_PI);
    return(sum);
  }
  double calc_nabla1(const double q) {
    // Computes the nabla xi_L integral, with a smoothing, since this
    // is numerically tricky.
    const double Rsmth2=1.0*1.0;
    const int Nk=kLin.size();
    const int Nint=10000;
    const double hh=(kLin[Nk-1]-kLin[0])/Nint;
    double sum=0;
#pragma omp parallel for reduction(+:sum)
    for (int i=1; i<Nint; ++i) {
      double xx = kLin[0]+i*hh;
      double kk = exp(xx);
      double k2 = kk*kk;
      int    jj = (int)(i*hh*dkinv);
      if (jj>=pLin.size()-2) jj=pLin.size()-2;
      double pk = exp(pLin[jj]+(xx-kLin[jj])*
                     (pLin[jj+1]-pLin[jj])/(kLin[jj+1]-kLin[jj]));
      double wk = exp(-k2*Rsmth2);
      std::vector<double> jl=sphBess(kk*q);
      int wt= 2+2*(i%2);
      sum  +=-kk*pk*(k2*k2*q*jl[1])*wt*wk;	// xi_lin'
    }
    sum *= hh/3.0/(2*M_PI*M_PI);
    return(sum);
  }
  double calc_nabla2(const double q) {
    // Computes the nabla^2 xi_L integral, with a smoothing, since this
    // is numerically tricky.
    const double Rsmth2=1.0*1.0;
    const int Nk=kLin.size();
    const int Nint=10000;
    const double hh=(kLin[Nk-1]-kLin[0])/Nint;
    double sum=0;
#pragma omp parallel for reduction(+:sum)
    for (int i=1; i<Nint; ++i) {
      double xx = kLin[0]+i*hh;
      double kk = exp(xx);
      double k2 = kk*kk;
      int    jj = (int)(i*hh*dkinv);
      if (jj>=pLin.size()-2) jj=pLin.size()-2;
      double pk = exp(pLin[jj]+(xx-kLin[jj])*
                     (pLin[jj+1]-pLin[jj])/(kLin[jj+1]-kLin[jj]));
      double wk = k2*exp(-k2*Rsmth2);
      std::vector<double> jl=sphBess(kk*q);
      int wt= 2+2*(i%2);
      sum  += kk*pk*(k2*jl[0])*wt*wk;		// \nabla^2 xi_lin
    }
    sum *= hh/3.0/(2*M_PI*M_PI);
    return(sum);
  }
  std::vector<double> calc_Jn(const double q) {
    // Computes the \mathcal{J}_n integrals, which are used in the shear
    // terms, and the other shear-related terms.
    const int Nk=kLin.size();
    const int Nint=10000;
    const double hh=(kLin[Nk-1]-kLin[0])/Nint;
    double sum1=0,sum2=0,sum3=0,sum4=0,sum5=0,sum6=0,sum7=0,sum8=0,sum9=0;
#pragma omp parallel for reduction(+:sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9)
    for (int i=1; i<Nint; ++i) {
      double xx = kLin[0]+i*hh;
      double ap = cos(M_PI/2.*exp(xx-kLin[Nk-1]));
      double kk = exp(xx);
      double k2 = kk*kk;
      double kq = kk*q;
      int    jj = (int)(i*hh*dkinv);
      if (jj>=pLin.size()-2) jj=pLin.size()-2;
      double pk = exp(pLin[jj]+(xx-kLin[jj])*
                     (pLin[jj+1]-pLin[jj])/(kLin[jj+1]-kLin[jj]));
      std::vector<double> jl=sphBess(kq);
      double j0,j1,j2,j3,j4;
      j0=jl[0];
      if (kq<0.9) {
        double kq2 = kq*kq;
        j1 = kq *(1./3.+ kq2*(-1./30.+kq2*(1./840.-kq2/45360.)));
        j2 = kq2*(1./15.+kq2*(-1./210.+kq2*(1./7560-kq2/498960.)));
        j3 = kq *kq2*(1./105.+kq2*(-1./1890+kq2*(1./83160-kq2/6486480.)));
        j4 = kq2*kq2*(1./945.+kq2*(-1./20790.+kq2/1081080.));
      }
      else {
        j1 =    jl[1]*kq;
        j2 = 3.*jl[1]  -jl[0];
        j3 = 5.*j2/(kq)-j1;
        j4 = 7.*j3/(kq)-j2;
      }
      int wt= 2+2*(i%2);
      sum1 += k2*pk*kk*(j2)*wt;
      sum2 += k2*pk*(2./15.*j1-1./5.*j3)*wt * ap;
      sum3 += k2*pk*(-1./5.*j1-1./5.*j3)*wt;
      sum4 += k2*pk*(j3)*wt;
      sum5 += k2*pk*kk*(-14*j0-40*j2+9*j4)/315.*wt*ap;
      sum6 += k2*pk*kk*(  7*j0+10*j2+3*j4)/105.*wt*ap;
      sum7 += k2*pk*kk*(        4*j2-3*j4)/ 21.*wt*ap;
      sum8 += k2*pk*kk*(       -3*j2-3*j4)/ 21.*wt*ap;
      sum9 += k2*pk*kk*(               j4)     *wt*ap;
    }
    sum5 *= hh/3.0/(2*M_PI*M_PI);
    sum6 *= hh/3.0/(2*M_PI*M_PI);
    sum7 *= hh/3.0/(2*M_PI*M_PI);
    sum8 *= hh/3.0/(2*M_PI*M_PI);
    sum9 *= hh/3.0/(2*M_PI*M_PI);
    double zeta= sum5*( 9*sum5+12*sum6+12*sum7+ 8*sum8+ 2*sum9)+
                 sum6*(        24*sum6+ 8*sum7+32*sum8+ 4*sum9)+
                 sum7*(               + 8*sum7+16*sum8+ 4*sum9)+
                 sum8*(                        24*sum8+ 8*sum9)+
                 sum9*(                                   sum9);
    std::vector<double> sum(8);
    sum[1] = sum1 * hh/3.0/(2*M_PI*M_PI);	// mathcal{J}_1
    sum[2] = sum2 * hh/3.0/(2*M_PI*M_PI);	// mathcal{J}_2
    sum[3] = sum3 * hh/3.0/(2*M_PI*M_PI);	// mathcal{J}_3
    sum[4] = sum4 * hh/3.0/(2*M_PI*M_PI);	// mathcal{J}_4
    sum[5] = 4    * sum[1]*sum[2];		// V_i^{12}
    sum[6] = 4./3.* sum[1]*sum[1];		// chi12
    sum[7] = 2*zeta;				// zeta
    return(sum);
  }
  void tabulateQfuncs() {
    // Finally, tabulate sigma2, etaPer, etaPar, etc.
    qmin = 1.0/exp(kLin[NkTable-1]);  if (qmin<0.2) qmin=0.2;
    qmax = 1.0/exp(kLin[    0    ]);  if (qmax>250) qmax=250;
    // First compute it on a coarse grid.
    std::vector<double> qvals;
    const int Nsample=150;
    try {
       qvals.resize(Nsample);
      etaPer.resize(Nsample);
      etaPar.resize(Nsample);
        uVal.resize(Nsample);
       xiLin.resize(Nsample);
      dxiLin.resize(Nsample);
       n2Lin.resize(Nsample);
       J2Lin.resize(Nsample);
       J3Lin.resize(Nsample);
       J4Lin.resize(Nsample);
      V12Lin.resize(Nsample);
      chiLin.resize(Nsample);
      zetLin.resize(Nsample);
    } catch(std::exception& e) {myexception(e);}
    delta=(qmax-qmin)/(Nsample-1);
    for (int i=0; i<Nsample; ++i) {
      double qq = qmin+i*delta;
      std::vector<double> qf=calcQfuncs(qq);
      std::vector<double> Jn=calc_Jn(qq);
       qvals[i] = qq;
      etaPer[i] = qf[0];
      etaPar[i] = qf[1];
      uVal[  i] = qf[2];
       xiLin[i] = qf[3] * qq*qq;
      dxiLin[i] = calc_nabla1(qq);
       n2Lin[i] = calc_nabla2(qq);
       J2Lin[i] = Jn[2];
       J3Lin[i] = Jn[3];
       J4Lin[i] = Jn[4];
      V12Lin[i] = Jn[5];
      chiLin[i] = Jn[6];
      zetLin[i] = Jn[7];
    }
    // then fit splines and retabulate it onto a finer grid.
    Spline etaPerSpline(qvals,etaPer);
    Spline etaParSpline(qvals,etaPar);
    Spline   uValSpline(qvals,uVal);
    Spline  xiLinSpline(qvals, xiLin);
    Spline dxiLinSpline(qvals,dxiLin);
    Spline  n2LinSpline(qvals, n2Lin);
    Spline  J2LinSpline(qvals, J2Lin);
    Spline  J3LinSpline(qvals, J3Lin);
    Spline  J4LinSpline(qvals, J4Lin);
    Spline V12LinSpline(qvals,V12Lin);
    Spline chiLinSpline(qvals,chiLin);
    Spline zetLinSpline(qvals,zetLin);
    try {
      etaPer.resize(NqTable);
      etaPar.resize(NqTable);
        uVal.resize(NqTable);
       xiLin.resize(NqTable);
      dxiLin.resize(NqTable);
       n2Lin.resize(NqTable);
       J2Lin.resize(NqTable);
       J3Lin.resize(NqTable);
       J4Lin.resize(NqTable);
      V12Lin.resize(NqTable);
      chiLin.resize(NqTable);
      zetLin.resize(NqTable);
    } catch(std::exception& e) {myexception(e);}
    sigma2 = calcSigma2();
    delta=(qmax-qmin)/(NqTable-1);
    for (int i=0; i<NqTable; ++i) {
      double qq = qmin+i*delta;
      etaPer[i] = etaPerSpline(qq);
      etaPar[i] = etaParSpline(qq);
      uVal[  i] =   uValSpline(qq);
       xiLin[i] =  xiLinSpline(qq)/qq/qq;
      dxiLin[i] = dxiLinSpline(qq);
       n2Lin[i] =  n2LinSpline(qq);
       J2Lin[i] =  J2LinSpline(qq);
       J3Lin[i] =  J3LinSpline(qq);
       J4Lin[i] =  J4LinSpline(qq);
      V12Lin[i] = V12LinSpline(qq);
      chiLin[i] = chiLinSpline(qq);
      zetLin[i] = zetLinSpline(qq);
    }
  }
  std::vector<double> interpQfuncs(const double q) {
    // Does a linear interpolation to return etaPer and etaPar.
    // Also returns U(q) and xi_lin, and some of the 1D integrals
    // needed for the shear terms.
    std::vector<double> qf(12);
    int k=(NqTable-1)*(q-qmin)/(qmax-qmin);
    if (q>qmin && q<qmax) {
      double dq = (q-(qmin+k*delta))/delta;
      qf[ 0]=etaPer[k]+dq*(etaPer[k+1]-etaPer[k]);
      qf[ 1]=etaPar[k]+dq*(etaPar[k+1]-etaPar[k]);
      qf[ 2]=  uVal[k]+dq*(  uVal[k+1]-  uVal[k]);
      qf[ 3]= xiLin[k]+dq*( xiLin[k+1]- xiLin[k]);
      qf[ 4]=dxiLin[k]+dq*(dxiLin[k+1]-dxiLin[k]);
      qf[ 5]= n2Lin[k]+dq*( n2Lin[k+1]- n2Lin[k]);
      qf[ 6]= J2Lin[k]+dq*( J2Lin[k+1]- J2Lin[k]);
      qf[ 7]= J3Lin[k]+dq*( J3Lin[k+1]- J3Lin[k]);
      qf[ 8]= J4Lin[k]+dq*( J4Lin[k+1]- J4Lin[k]);
      qf[ 9]=V12Lin[k]+dq*(V12Lin[k+1]-V12Lin[k]);
      qf[10]=chiLin[k]+dq*(chiLin[k+1]-chiLin[k]);
      qf[11]=zetLin[k]+dq*(zetLin[k+1]-zetLin[k]);
    }
    else {
      const double TINY=1e-10;
      if (q<qmin) {
        qf[ 0]=sigma2 - TINY;
        qf[ 1]=sigma2 - TINY;
        qf[ 2]=0;
        qf[ 3]=0;
        qf[ 4]=0;
        qf[ 5]=0;
        qf[ 6]=0;
        qf[ 7]=0;
        qf[ 8]=0;
        qf[ 9]=0;
        qf[10]=0;
        qf[11]=0;
      }
      if (q>qmax) {
        qf[ 0]=etaPer[NqTable-1];
        qf[ 1]=etaPar[NqTable-1];
        qf[ 2]=  uVal[NqTable-1];
        qf[ 3]= xiLin[NqTable-1];
        qf[ 4]=dxiLin[NqTable-1];
        qf[ 5]= n2Lin[NqTable-1];
        qf[ 6]= J2Lin[NqTable-1];
        qf[ 7]= J3Lin[NqTable-1];
        qf[ 8]= J4Lin[NqTable-1];
        qf[ 9]=V12Lin[NqTable-1];
        qf[10]=chiLin[NqTable-1];
        qf[11]=zetLin[NqTable-1];
      }
    }
    return(qf);
  }
  std::vector<double> calcAmat(const double q[]) {
    // Returns the 3x3 matrix A (Eq. 28 of CLPT).
    double qhat[3],qq=0;
    for (int i=0; i<3; ++i) qq += q[i]*q[i];  qq=sqrt(qq);
    for (int i=0; i<3; ++i) qhat[i] = q[i]/qq;
    std::vector<double> qf=interpQfuncs(qq);
    double F =2*(sigma2-qf[0]);	// sigma_perp^2
    double G =2*(qf[0] -qf[1]);	// sigma_par ^2 - sigma_perp^2
    std::vector<double> Amat(9);
    for (int i=0; i<3; ++i)
      for (int j=i; j<3; ++j) {
        Amat[3*i+j] = Amat[3*j+i]  = G*qhat[i]*qhat[j];
        if (i==j)     Amat[3*i+i] += F;
      }
    return(Amat);
  }
  std::vector<double> calcAinv(const double q[]) {
    // Returns the inverse of the 3x3 matrix A (Eq. 28 of CLPT).
    // Also returns the determinant of Ainv as the last (extra) element.
    // The Sherman-Morrison formula says that
    //   (A+b.cT)^{-1}=Ainv - Ainv.b.cT.Ainv/[1+cT.Ainv.b]
    // so our (F.1+G.q.q)^{-1}=1/F-G.q.q/F/[F+G]
    // For real-space (r-q).Ainv.(r-q) depends only on r, q and rq.mu.
    // Moving into redshift space simply requires us to
    // divide the zz element by (1+f)^2 [and change det], but now we
    // need to do the phi integral numerically as well.
    double qhat[3],qq=0;
    for (int i=0; i<3; ++i) qq += q[i]*q[i];  qq=sqrt(qq);
    for (int i=0; i<3; ++i) qhat[i] = q[i]/qq;
    std::vector<double> qf=interpQfuncs(qq);
    double  F =2*(sigma2-qf[0]);	// sigma_perp^2
    double  G =2*(qf[0] -qf[1]);	// sigma_par ^2 - sigma_perp^2
    double FpG=2*(sigma2-qf[1]);	// sigma_par ^2
    std::vector<double> Ainv(10);
    for (int i=0; i<3; ++i)
      for (int j=i; j<3; ++j) {
        Ainv[3*i+j] = Ainv[3*j+i]  = -G*qhat[i]*qhat[j]/F/FpG;
        if (i==j)     Ainv[3*i+i] += 1.0/F;
      }
    // Now set detA.  Use det(cM)=c^n det(M) and det(I+u.vT)=1+uT.v so that
    // det(F.delta_ij+G.qhat_i.qhat_j) = F^3(1+G/F) = F^2(F+G).
    // Also note that F+G is 2[sigma^2-etaPar] for our case, which is the
    // same thing as sigma_{||}, while F is \sigma_\perp, thus
    // detA=(sigma_perp^2.sigma_par)^2
    Ainv[9] = 1.0/(F*F*FpG);
    return(Ainv);
  }
  double zeldovichIntegrand(const double r[], const double q[], const double f){
    // Returns the integrand for the 1+\xi integral (Eq. 34 of CLPT).
    const double twoPi3=248.05021344239853;
    std::vector<double> Ainv=calcAinv(q);	// Also returns detAinv.
    if (f>0) { // If we are in redshift space, need to correct for A->RAR.
      for (int i=0; i<3; ++i) {
        Ainv[3*i+2] /= (1+f);
        Ainv[3*2+i] /= (1+f);
      }
      Ainv[9] /= (1+f)*(1+f);
    }
    double res=0;
    for (int i=0; i<3; ++i)
      for (int j=0; j<3; ++j)
        res += (r[i]-q[i])*Ainv[3*i+j]*(r[j]-q[j]);
    res = exp(-0.5*res)*sqrt(Ainv[9]/twoPi3);
    return(res);
  }
public:
  void init(const char fname[]) {
    // Initialize the G-L integration points and weights.
    gl.set(128);
    // Load the linear power spectrum from a file, and expand it out
    // if necessary.
    readPowerSpectrum(fname);
    // Finally pre-tabulate the q-dependent functions.
    tabulateQfuncs();
    // and the normalization term for EFT.
    eftNorm = calcEFTnorm();
  } // End init.
  Zeldovich() {}
  Zeldovich(const char fname[]) {
    init(fname);
  }
  ~Zeldovich() {}
  void print_eta() {
    // A convenience/debug feature to print sigma^2, etaPerp and etaPar.
    // Also prints useful combinations of these, and prints U(q).
    std::cout<<"sigma^2="<<std::scientific<<sigma2<<std::endl;
    std::cout<<"# "<<std::setw(10)<<"q(Mpc/h)"
                   <<std::setw(12)<<"Eta_perp"
                   <<std::setw(12)<<"Eta_par"
                   <<std::setw(12)<<"Sig2_perp"
                   <<std::setw(12)<<"Sig2_par"
                   <<std::setw(12)<<"U"
                   << std::endl;
    for (int i=0; i<NqTable; ++i) {
      double qq     = qmin + i*delta;
      double sigper = 2*(sigma2-etaPer[i]);
      double sigpar = 2*(sigma2-etaPar[i]);
      std::cout
        <<std::scientific<<std::setw(12)<<std::setprecision(4)<<qq
        <<std::scientific<<std::setw(12)<<std::setprecision(4)<<etaPer[i]
        <<std::scientific<<std::setw(12)<<std::setprecision(4)<<etaPar[i]
        <<std::scientific<<std::setw(12)<<std::setprecision(4)<<sigper
        <<std::scientific<<std::setw(12)<<std::setprecision(4)<<sigpar
        <<std::scientific<<std::setw(12)<<std::setprecision(4)<<uVal[i]
        <<std::endl;
    }
    std::cout.flush();
  }
  double xiL(const double r) {
    // The real-space, linear correlation function at r.
    // This is not tested for very large or small values of r.
    std::vector<double> qf=interpQfuncs(r);
    return(qf[3]);
  }
  double xiZ(const double rval) {
    // Returns the real-space, matter, Zel'dovich correlation function at r.
    // This is not tested for very large or small values of r.
    // The integration is over x=q-r, in length and angle with the
    // azimuthal integral being trivial.
    const double xmin=0;
    const double xmax=10*sqrt(sigma2);
    const double rr[3]={0,0,rval};
    const double r2   =rval*rval;
    const int    Nx=2000;
    const double dx=(xmax-xmin)/Nx;
    double xi=0;
    for (int ixx=0; ixx<Nx; ++ixx) {
      double xx=xmin+(ixx+0.5)*dx;
      double x2=xx*xx;
      for (int imu=0; imu<gl.N; ++imu) {
        double mu = gl.x[imu];
        // Compute vec{q}=vec{r}+vec{x} with vec{r}=r.zhat,
        // so r_z=r, x_z=x*mu, cos_rq=(r_z+x_z)/q.
        double qlen = sqrt(r2+x2+2*rval*xx*mu);
        double qcos = (rval+xx*mu)/qlen;
        double qq[3] ={qlen*sqrt(1-qcos*qcos),0,qlen*qcos};
        if (qlen>qmin && qlen<qmax)
          xi += x2 * zeldovichIntegrand(rr,qq,0) * gl.w[imu];
      }
    }
    xi *= dx;		// Convert sum to integral.
    xi *= 2*M_PI;	// The azimuthal integral.
    xi -= 1.0;		// Calculated 1+xi, subtract 1 for xi.
    return(xi);
  }
  std::vector<double> xiContributions(const double rval, const double Aeft) {
    // Returns the different contributions to the real-space Zel'dovich
    // correlation function for locally biased tracers.
    // This is not tested for very large or small values of r.
    // The integration is over x=q-r, in length and angle with the
    // azimuthal integral being trivial.
    const double xmin=0;
    const double xmax=10*sqrt(sigma2);
    const double rr[3]={0,0,rval};
    const double r2   =rval*rval;
    const int    Nx=500;
    const double dx=(xmax-xmin)/Nx;
    std::vector<double> xi(6);
    for (int ixx=0; ixx<Nx; ++ixx) {
      double xx=xmin+(ixx+0.5)*dx;
      double x2=xx*xx;
      for (int imu=0; imu<gl.N; ++imu) {
        double mu = gl.x[imu];
        // Compute vec{q}=vec{r}+vec{x} with vec{r}=r.zhat,
        // so r_z=r, x_z=x*mu, cos_rq=(r_z+x_z)/q.
        double qlen = sqrt(r2+x2+2*rval*xx*mu);
        double qcos = (rval+xx*mu)/qlen;
        double qsin = sqrt(1-qcos*qcos);
        double qq[3] ={qlen*qsin,0,qlen*qcos};
        double qh[3] ={     qsin,0,     qcos};
        if (qlen>qmin && qlen<qmax) {
          // For the unbiased tracers we only need this--all other terms
          // are multiplied by this anyway.
          double pref = x2 * zeldovichIntegrand(rr,qq,0) * gl.w[imu];
          // For the bias terms, compute U,xi and Ainv (even though in above).
          std::vector<double> qf  =interpQfuncs(qlen);
          std::vector<double> Ainv=calcAinv(qq);
          // Construct the auxilliary matrix/vectors g, G of CLPT Eq. (45)
          double g[3],G[9];
          for (int i=0; i<3; ++i) {
            g[i]=0;
            for (int j=0; j<3; ++j)
              g[i] += Ainv[3*i+j]*(qq[j]-rr[j]);
          }
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              G[3*i+j]=Ainv[3*i+j]-g[i]*g[j];
          double Ug,UUG,trG;  Ug=UUG=trG=0;
          for (int i=0; i<3; ++i) {
            Ug += (qf[2]*qh[i])*g[i];
            trG+= G[3*i+i];
            for (int j=0; j<3; ++j)
              UUG += qf[2]*qf[2]*qh[i]*qh[j]*G[3*i+j];
          }
          // Now do the 1, Fp, Fpp, Fp^2, Fp.Fpp, Fpp^2 terms.
          xi[0] +=    pref*(1+Aeft*trG*eftNorm);
          xi[1] += -2*pref*Ug;
          xi[2] +=   -pref*UUG;
          xi[3] +=    pref*(qf[3]-UUG);
          xi[4] += -2*pref*qf[3]*Ug;
          xi[5] +=0.5*pref*qf[3]*qf[3];
        }
      }
    }
    for (int j=0; j<xi.size(); ++j) {
      xi[j] *= dx;	// Convert sum to integral.
      xi[j] *= 2*M_PI;	// The azimuthal integral.
    }
    xi[0] -= 1.0;	// Calculated 1+xi, subtract 1 for xi.
    return(xi);
  }
  std::vector<double> v12(const double rval) {
    // Returns the different contributions to the mean infall velocity
    // for locally biased tracers.  Only the line-of-sight component is
    // returned and the result should be multiplied by f and divided by
    // 1+xi(real).
    // This is not tested for very large or small values of r.
    // The integration is over x=q-r, in length and angle with the
    // azimuthal integral being trivial.
    const double xmin=0;
    const double xmax=10*sqrt(sigma2);
    const double rr[3]={0,0,rval};
    const double r2   =rval*rval;
    const int    Nx=500;
    const double dx=(xmax-xmin)/Nx;
    std::vector<double> vv(6);
    for (int ixx=0; ixx<Nx; ++ixx) {
      double xx=xmin+(ixx+0.5)*dx;
      double x2=xx*xx;
      for (int imu=0; imu<gl.N; ++imu) {
        double mu = gl.x[imu];
        // Compute vec{q}=vec{r}+vec{x} with vec{r}=r.zhat,
        // so r_z=r, x_z=x*mu, cos_rq=(r_z+x_z)/q.
        double qlen = sqrt(r2+x2+2*rval*xx*mu);
        double qcos = (rval+xx*mu)/qlen;
        double qsin = sqrt(1-qcos*qcos);
        double qq[3] ={qlen*qsin,0,qlen*qcos};
        double qh[3] ={     qsin,0,     qcos};
        if (qlen>qmin && qlen<qmax) {
          // For the unbiased tracers we only need this--all other terms
          // are multiplied by this anyway.
          double pref = x2 * zeldovichIntegrand(rr,qq,0) * gl.w[imu];
          // For the bias terms compute U, xi and Ainv (even though in above).
          std::vector<double> qf  =interpQfuncs(qlen);
          std::vector<double> Ainv=calcAinv(qq);
          std::vector<double> Amat=calcAmat(qq);
          // Construct the auxilliary matrix/vectors g, G of CLPT Eq. (45)
          double g[3],U[3],G[9];
          for (int i=0; i<3; ++i) {
            g[i]=0;
            for (int j=0; j<3; ++j)
              g[i] += Ainv[3*i+j]*(qq[j]-rr[j]);
            U[i]=qf[2]*qh[i];
          }
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              G[3*i+j]=Ainv[3*i+j]-g[i]*g[j];
          double Ug,gA,UGA;  Ug=gA=UGA=0;
          for (int i=0; i<3; ++i) {
            gA += g[i]*Amat[3*2+i];
            Ug += U[i]*g[i];
            for (int j=0; j<3; ++j)
              UGA += U[i]*G[3*i+j]*Amat[3*2+j];
          }
          // Now do the 1, Fp, Fpp, Fp^2, Fp.Fpp, Fpp^2 terms.
          vv[0] +=   -pref * gA;
          vv[1] +=  2*pref *(U[2]-UGA);
          vv[2] += -2*pref *Ug*U[2];
          vv[3] +=   -pref *(2*Ug*U[2]+qf[3]*gA);
          vv[4] +=  2*pref *qf[3]*U[2];
          vv[5] +=  0;
        }
      }
    }
    for (int j=0; j<vv.size(); ++j) {
      vv[j] *= dx;	// Convert sum to integral.
      vv[j] *= 2*M_PI;	// The azimuthal integral.
    }
    return(vv);
  }
  std::vector<double> s12(const double rval) {
    // Returns the different contributions to the velocity dispersion
    // for locally biased tracers.  Both sigma_perp and sigma_par are
    // returned and the result should be multiplied by f^2 and divided by
    // 1+xi(real).  Beware the ordering here!!
    // This is not tested for very large or small values of r.
    // The integration is over x=q-r, in length and angle with the
    // azimuthal integral being trivial.
    const double xmin=0;
    const double xmax=10*sqrt(sigma2);
    const double rr[3]={0,0,rval};
    const double r2   =rval*rval;
    const int    Nx=500;
    const double dx=(xmax-xmin)/Nx;
    std::vector<double> ss(12);
    for (int ixx=0; ixx<Nx; ++ixx) {
      double xx=xmin+(ixx+0.5)*dx;
      double x2=xx*xx;
      for (int imu=0; imu<gl.N; ++imu) {
        double mu = gl.x[imu];
        // Compute vec{q}=vec{r}+vec{x} with vec{r}=r.zhat,
        // so r_z=r, x_z=x*mu, cos_rq=(r_z+x_z)/q.
        double qlen = sqrt(r2+x2+2*rval*xx*mu);
        double qcos = (rval+xx*mu)/qlen;
        double qsin = sqrt(1-qcos*qcos);
        double qq[3] ={qlen*qsin,0,qlen*qcos};
        double qh[3] ={     qsin,0,     qcos};
        if (qlen>qmin && qlen<qmax) {
          // For the unbiased tracers we only need this--all other terms
          // are multiplied by this anyway.
          double pref = x2 * zeldovichIntegrand(rr,qq,0) * gl.w[imu];
          // For the bias terms compute U, xi and Ainv (even though in above).
          std::vector<double> qf  =interpQfuncs(qlen);
          std::vector<double> Ainv=calcAinv(qq);
          std::vector<double> Amat=calcAmat(qq);
          // Construct the auxilliary matrix/vectors g, G of CLPT Eq. (45)
          double g[3],U[3],G[9];
          for (int i=0; i<3; ++i) {
            g[i]=0;
            for (int j=0; j<3; ++j)
              g[i] += Ainv[3*i+j]*(qq[j]-rr[j]);
            U[i]=qf[2]*qh[i];
          }
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              G[3*i+j]=Ainv[3*i+j]-g[i]*g[j];
          // Work out these dot products for n=m=2, i.e. along z, for sigma_par.
          double Ug,U2,gA,gAU,AGA,trA;  Ug=U2=gA=gAU=AGA=trA=0;
          for (int i=0; i<3; ++i) {
            gA += g[i]*Amat[3*2+i];
            Ug += U[i]*g[i];
            U2 += U[i]*U[i];
            trA+= Amat[3*i+i];
            for (int j=0; j<3; ++j) {
              gAU += g[i]*Amat[3*i+j]*U[j];
              AGA += Amat[3*2+i]*G[3*i+j]*Amat[3*2+j];
            }
          }
          // Now do the 1, Fp, Fpp, Fp^2, Fp.Fpp, Fpp^2 terms for \sigma_par^2.
          ss[ 0] +=    pref *(Amat[3*2+2]-AGA);
          ss[ 1] += -2*pref *(2*gA*U[2]+Ug*Amat[3*2+2]);
          ss[ 2] +=  2*pref * U[2]*U[2];
          ss[ 3] +=    pref *(2*U[2]*U[2]+qf[3]*Amat[3*2+2]);
          ss[ 4]  =  0;
          ss[ 5]  =  0;
          // Next work out the trace components, i.e. summed over n=m.
          AGA=0;
          for (int m=0; m<3; ++m) {
            for (int i=0; i<3; ++i) {
              for (int j=0; j<3; ++j) {
                AGA += Amat[3*m+i]*G[3*i+j]*Amat[3*j+m];
              }
            }
          }
          ss[ 6] +=    pref *(trA-AGA);
          ss[ 7] += -2*pref *(2*gAU+Ug*trA);
          ss[ 8] +=  2*pref * U2;
          ss[ 9] +=    pref *(2*U2+qf[3]*trA);
          ss[10]  =  0;
          ss[11]  =  0;
        }
      }
    }
    // Now sigma_perp is related to the trace by s_perp^2=(1/2)[Tr-sig_par^2]
    for (int j=6; j<ss.size(); ++j)
      ss[j] = 0.5*(ss[j] - ss[j-6]);
    for (int j=0; j<ss.size(); ++j) {
      ss[j] *= dx;	// Convert sum to integral.
      ss[j] *= 2*M_PI;	// The azimuthal integral.
    }
    // Note we return parallel then perpendicular/transverse!
    return(ss);
  }
};





class LSM: public Zeldovich {
// Implements the Lagrangian streaming model.
// Currently set up to take a linear power spectrum, f, F1 and F2
// at input, but could be modified to allow f, F1 and F2 to vary.
protected:
  Spline		xispl,vvspl,stspl,spspl;
  Spline		R1spl,R2spl,Q1spl,Q2spl,Q5spl,Q8spl,Qsspl;
  std::vector<double>	X11,X22,X13,Y11,Y22,Y13,X1210,Y1210;
  std::vector<double>	V1,V3,TT;
  std::vector<double>	U1,U3,U220,U211,S2D;
  LPT			lpt;
  std::vector<double> calcEfuncs(const double q) {
    // Computes the "extra" functions of Q which come in beyond the
    // Zel'dovich approximation.
    const int Nk=kLin.size();
    const double kmax=exp(kLin[Nk-1]);
    int Nint=(int)(8*kmax*q+512);
    if (Nint>=20000) Nint=20000;
    const double hh=(kLin[Nk-1]-kLin[0])/Nint;
    double sum0=0,sum1=0,sum2=0,sum3=0,sum4=0,sum5=0,sum6=0,sum7=0,
           sum8=0,sum9=0,sum10=0,sum11=0,sum12=0,sumS=0;
#pragma omp parallel for reduction(+:sum0,sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9,sum10,sum11,sum12,sumS)
    for (int i=1; i<Nint; ++i) {
      double xx = kLin[0]+i*hh;
      double kk = exp(xx);
      double k2 = kk*kk;
      double kq = kk*q;
      double R1 = R1spl(kk);
      double R2 = R2spl(kk);
      double Q1 = Q1spl(kk);
      double Q2 = Q2spl(kk);
      double Q5 = Q5spl(kk);
      double Q8 = Q8spl(kk);
      double Qs = Qsspl(kk);
      double ap = cos(M_PI/2.*kk/kmax);
      std::vector<double> jl=sphBess(kk*q);
      double j1 = kq*jl[1];
      double j2,j3;
      if (kq<0.1) {
        j2 = pow(kq,2.0)/15.  - pow(kq,4.0)/210.;
        j3 = pow(kq,3.0)/105. - pow(kq,5.0)/1890.;
      }
      else {
        j2 = 3.*jl[1]-jl[0];
        j3 = 5.*j2/(kq)-(kq)*jl[1];
      }
      int wt= 2+2*(i%2);
      sum0 += kk*wt*(9./98.*Q1*(2./3.-2*jl[1]));		// X^{(22)}
      sum1 += kk*wt*(5./21.*R1*(2./3.-2*jl[1]));		// X^{(13)}
      sum2 += kk*wt*(9./98.*Q1*(-2*jl[0]+6*jl[1]));		// Y^{(22)}
      sum3 += kk*wt*(5./21.*R1*(-2*jl[0]+6*jl[1]));		// Y^{(13)}
      sum4 += kk*wt*(2*(R1-R2)+3*R1*jl[0]-3*(3*R1+4*R2+2*Q5)*jl[1])/14.;//X1210
      sum5 += kk*wt*(3*R1+4*R2+2*Q5)*(jl[0]-3*jl[1])*(-3./14.);// Y_{10}^{(12)}
      sum6 +=    wt*(R1*j1)*(-3./7.);			// V_1^{(112)}
      sum7 +=    wt*(Q1*j1)*(-3./7.);			// V_3^{(112)}
      sumS +=    wt*(2*R1+4*R2+Q1+2*Q2)*(3./7.*j2/(kk*q));	// S^{(112)}
      sum8 +=    wt*(2*R1+4*R2+Q1+2*Q2)*j3*(-3./7.);	// T^{(112)}
      sum9 += k2*wt*(R1*j1)*(-5./21.);			// U^{(3)}
      sum10+= k2*wt*(Q8*j1)*(-3./7.);			// U_{20}^{(2)}
      sum11+= k2*wt*((R1+R2)*j1)*(-6./7.);		// U_{11}^{(2)}
      sum12+= k2*wt*(Qs*j1)*(-2./7.)*ap;		// Shear term
    }
    sum6 += sumS;
    sum7 += sumS;
    std::vector<double> sum(16);
    sum[ 1] = sum0 * hh/3.0/(2*M_PI*M_PI);
    sum[ 2] = sum1 * hh/3.0/(2*M_PI*M_PI);
    sum[ 4] = sum2 * hh/3.0/(2*M_PI*M_PI);
    sum[ 5] = sum3 * hh/3.0/(2*M_PI*M_PI);
    sum[ 6] = sum4 * hh/3.0/(2*M_PI*M_PI);
    sum[ 7] = sum5 * hh/3.0/(2*M_PI*M_PI);
    sum[ 8] = sum6 * hh/3.0/(2*M_PI*M_PI);
    sum[ 9] = sum7 * hh/3.0/(2*M_PI*M_PI);	
    sum[10] = sum8 * hh/3.0/(2*M_PI*M_PI);	
    sum[12] = sum9 * hh/3.0/(2*M_PI*M_PI);	
    sum[13] = sum10* hh/3.0/(2*M_PI*M_PI);	
    sum[14] = sum11* hh/3.0/(2*M_PI*M_PI);	
    sum[15] = sum12* hh/3.0/(2*M_PI*M_PI);	
    // Now tabulate the pieces going as Plin.
    sum0=sum1=sum2=0;
#pragma omp parallel for reduction(+:sum0,sum1,sum2)
    for (int i=1; i<Nint; ++i) {
      double xx = kLin[0]+i*hh;
      double kk = exp(xx);
      double k2 = kk*kk;
      double kq = kk*q;
      int    jj = (int)(i*hh*dkinv);
      if (jj>=pLin.size()-2) jj=pLin.size()-2;
      double pk = exp(pLin[jj]+(xx-kLin[jj])*
                     (pLin[jj+1]-pLin[jj])/(kLin[jj+1]-kLin[jj]));
      std::vector<double> jl=sphBess(kk*q);
      double j1 = kq*jl[1];
      int wt= 2+2*(i%2);
      sum0 += kk*wt*pk*(2./3.-2*jl[1]);		// X^{(11)}
      sum1 += kk*wt*pk*(-2.*jl[0]+6*jl[1]);	// Y^{(11)}
      sum2 += k2*wt*pk*(-j1);			// U^{(1)}
    }
    sum[ 0] = sum0 * hh/3.0/(2*M_PI*M_PI);
    sum[ 3] = sum1 * hh/3.0/(2*M_PI*M_PI);
    sum[11] = sum2 * hh/3.0/(2*M_PI*M_PI);
    return(sum);
  }
  void tabulateEfuncs() {
    // Tabulate the "extra" functions.
    // First compute them on a coarse grid.
    std::vector<double> qq;
    const int Nsample=150;
    try {
      qq.resize(Nsample);
      X11.resize(Nsample);
      X22.resize(Nsample);
      X13.resize(Nsample);
      Y11.resize(Nsample);
      Y22.resize(Nsample);
      Y13.resize(Nsample);
      X1210.resize(Nsample);
      Y1210.resize(Nsample);
      V1.resize(Nsample);
      V3.resize(Nsample);
      TT.resize(Nsample);
      U1.resize(Nsample);
      U3.resize(Nsample);
      U220.resize(Nsample);
      U211.resize(Nsample);
      S2D.resize(Nsample);
    } catch(std::exception& e) {myexception(e);}
    delta=(qmax-qmin)/(Nsample-1);
    for (int i=0; i<Nsample; ++i) {
      double q = qmin+i*delta;
      std::vector<double> ef=calcEfuncs(q);
      qq[i]   = q;
      X11[i]  = ef[ 0];
      X22[i]  = ef[ 1];
      X13[i]  = ef[ 2];
      Y11[i]  = ef[ 3];
      Y22[i]  = ef[ 4];
      Y13[i]  = ef[ 5];
      X1210[i]= ef[ 6];
      Y1210[i]= ef[ 7];
      V1[i]   = ef[ 8];
      V3[i]   = ef[ 9];
      TT[i]   = ef[10];
      U1[i]   = ef[11];
      U3[i]   = ef[12];
      U220[i] = ef[13];
      U211[i] = ef[14];
      S2D[i]  = ef[15];
    }
    // then fit splines and retabulate it onto a finer grid.
    Spline X11Spline(qq,X11);
    Spline X22Spline(qq,X22);
    Spline X13Spline(qq,X13);
    Spline Y11Spline(qq,Y11);
    Spline Y22Spline(qq,Y22);
    Spline Y13Spline(qq,Y13);
    Spline X1210Spline(qq,X1210);
    Spline Y1210Spline(qq,Y1210);
    Spline V1Spline(qq,V1);
    Spline V3Spline(qq,V3);
    Spline TTSpline(qq,TT);
    Spline U1Spline(qq,U1);
    Spline U3Spline(qq,U3);
    Spline U220Spline(qq,U220);
    Spline U211Spline(qq,U211);
    Spline S2DSpline(qq,S2D);
    try {
      X11.resize(NqTable);
      X22.resize(NqTable);
      X13.resize(NqTable);
      Y11.resize(NqTable);
      Y22.resize(NqTable);
      Y13.resize(NqTable);
      X1210.resize(NqTable);
      Y1210.resize(NqTable);
      V1.resize(NqTable);
      V3.resize(NqTable);
      TT.resize(NqTable);
      U1.resize(NqTable);
      U3.resize(NqTable);
      U220.resize(NqTable);
      U211.resize(NqTable);
      S2D.resize(NqTable);
    } catch(std::exception& e) {myexception(e);}
    delta=(qmax-qmin)/(NqTable-1);
    for (int i=0; i<NqTable; ++i) {
      double q  = qmin+i*delta;
      X11[i]  = X11Spline(q);
      X22[i]  = X22Spline(q);
      X13[i]  = X13Spline(q);
      Y11[i]  = Y11Spline(q);
      Y22[i]  = Y22Spline(q);
      Y13[i]  = Y13Spline(q);
      X1210[i]= X1210Spline(q);
      Y1210[i]= Y1210Spline(q);
      V1[i]   = V1Spline(q);
      V3[i]   = V3Spline(q);
      TT[i]   = TTSpline(q);
      U1[i]   = U1Spline(q);
      U3[i]   = U3Spline(q);
      U220[i] = U220Spline(q);
      U211[i] = U211Spline(q);
      S2D[i]  = S2DSpline(q);
    }
  }
  void writeSaveFile(const char fname[]) {
    // Save the XX, YY, etc. arrays to a file.
    std::ofstream fs(fname,std::ios::trunc);
    if (!fs) {
      std::cerr<<"Unable to open "<<fname<<" for writing."<<std::endl;
      myexit(1);
    }
    for (int i=0; i<NqTable; ++i)
      fs << std::scientific << std::setw(20) << std::setprecision(9) << X11[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << X22[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << X13[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << Y11[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << Y22[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << Y13[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << X1210[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << Y1210[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << V1[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << V3[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << TT[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << U1[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << U3[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << U220[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << U211[i]
         << std::scientific << std::setw(20) << std::setprecision(9) << S2D[i]
         << std::endl;
    fs.close();
  }
  void readSaveFile(const char fname[]) {
    // Read the XX, YY, etc. arrays from a file.
    try {
      X11.resize(NqTable);
      X22.resize(NqTable);
      X13.resize(NqTable);
      Y11.resize(NqTable);
      Y22.resize(NqTable);
      Y13.resize(NqTable);
      X1210.resize(NqTable);
      Y1210.resize(NqTable);
      V1.resize(NqTable);
      V3.resize(NqTable);
      TT.resize(NqTable);
      U1.resize(NqTable);
      U3.resize(NqTable);
      U220.resize(NqTable);
      U211.resize(NqTable);
      S2D.resize(NqTable);
    } catch(std::exception& e) {myexception(e);}
    std::ifstream fs(fname);
    if (!fs) {
      std::cerr<<"Unable to open "<<fname<<" for reading."<<std::endl;
      myexit(1);
    }
    for (int i=0; i<NqTable; ++i) {
      std::string ss;
      getline(fs,ss);
      if (fs.fail()) {
        std::cerr<<"Error reading line "<<i<<" of "<<fname<<std::endl;
      }
      std::istringstream(ss) >> X11[i] >> X22[i] >> X13[i]
                             >> Y11[i] >> Y22[i] >> Y13[i]
                             >> X1210[i] >> Y1210[i]
                             >> V1[i] >> V3[i] >> TT[i] >> U1[i] >> U3[i]
                             >> U220[i] >> U211[i] >> S2D[i];
    }
    fs.close();
  }
  std::vector<double> interpEfuncs(const double q) {
    // Does a linear interpolation to return the "extra" functions.
    std::vector<double> ef(16);
    int k=(NqTable-1)*(q-qmin)/(qmax-qmin);
    if (q>qmin && q<qmax) {
      double dq = (q-(qmin+k*delta))/delta;
      ef[ 0]=X11[k]+dq*(X11[k+1]-X11[k]);
      ef[ 1]=X22[k]+dq*(X22[k+1]-X22[k]);
      ef[ 2]=X13[k]+dq*(X13[k+1]-X13[k]);
      ef[ 3]=Y11[k]+dq*(Y11[k+1]-Y11[k]);
      ef[ 4]=Y22[k]+dq*(Y22[k+1]-Y22[k]);
      ef[ 5]=Y13[k]+dq*(Y13[k+1]-Y13[k]);
      ef[ 6]=X1210[k]+dq*(X1210[k+1]-X1210[k]);
      ef[ 7]=Y1210[k]+dq*(Y1210[k+1]-Y1210[k]);
      ef[ 8]=V1[k]+dq*(V1[k+1]-V1[k]);
      ef[ 9]=V3[k]+dq*(V3[k+1]-V3[k]);
      ef[10]=TT[k]+dq*(TT[k+1]-TT[k]);
      ef[11]=U1[k]+dq*(U1[k+1]-U1[k]);
      ef[12]=U3[k]+dq*(U3[k+1]-U3[k]);
      ef[13]=U220[k]+dq*(U220[k+1]-U220[k]);
      ef[14]=U211[k]+dq*(U211[k+1]-U211[k]);
      ef[15]=S2D[k]+dq*(S2D[k+1]-S2D[k]);
    }
    else {
      if (q<qmin) {
        ef[0]=ef[1]=ef[2]=ef[3]=ef[4]=ef[5]=ef[6]=ef[7]=ef[8]=ef[9]
             =ef[10]=ef[11]=ef[12]=ef[13]=ef[14]=ef[15]=0;
      }
      if (q>qmax) {
        ef[ 0]=X11[NqTable-1];
        ef[ 1]=X22[NqTable-1];
        ef[ 2]=X13[NqTable-1];
        ef[ 3]=Y11[NqTable-1];
        ef[ 4]=Y22[NqTable-1];
        ef[ 5]=Y13[NqTable-1];
        ef[ 6]=X1210[NqTable-1];
        ef[ 7]=Y1210[NqTable-1];
        ef[ 8]=V1[NqTable-1];
        ef[ 9]=V3[NqTable-1];
        ef[10]=TT[NqTable-1];
        ef[11]=U1[NqTable-1];
        ef[12]=U3[NqTable-1];
        ef[13]=U220[NqTable-1];
        ef[14]=U211[NqTable-1];
        ef[15]=S2D[NqTable-1];
      }
    }
    return(ef);
  }
  void old_setupQR() {
    // Set up the Q and R's that we need.  We zero pad these out
    // to the maximum kLin.
    const int NkTemp=250;
    std::vector<double> ka(NkTemp),R1(NkTemp),R2(NkTemp);
    std::vector<double> Q1(NkTemp),Q2(NkTemp),Q5(NkTemp),Q8(NkTemp);
    const double kswitch=6.0,kdamp=3.00;
    if (exp(kLin[kLin.size()-1])<=kswitch) {
      std::cout<<"Power spectrum file does not go to high enough k."<<std::endl;
      myexit(666);
    }
#pragma omp parallel for
    for (int i=1; i<NkTemp; ++i) {
      double kk = exp( kLin[0]+i*(log(kswitch)-kLin[0])/(NkTemp-1) );
      double dd = exp(-(kk/kdamp)*(kk/kdamp));
      std::vector<double> Qn=lpt.Qn(kk);
      std::vector<double> Rn=lpt.Rn(kk);
      ka[i]=kk; R1[i]=Rn[0]*dd; R2[i]=Rn[1]*dd;
      Q1[i]=Qn[1]*dd; Q2[i]=Qn[2]*dd; Q5[i]=Qn[5]*dd; Q8[i]=Qn[8]*dd;
    }
    if (exp(kLin[kLin.size()-1])>kswitch)
      for (int i=0; i<8; ++i) {
        double kk = kswitch + (i+1)*(exp(kLin[kLin.size()-1])-kswitch)/8.0;
        ka.push_back(kk);
        R1.push_back(0.);
        R2.push_back(0.);
        Q1.push_back(0.);
        Q2.push_back(0.);
        Q5.push_back(0.);
        Q8.push_back(0.);
      }
    // and fit splines to them.
    R1spl.init(ka,R1);
    R2spl.init(ka,R2);
    Q1spl.init(ka,Q1);
    Q2spl.init(ka,Q2);
    Q5spl.init(ka,Q5);
    Q8spl.init(ka,Q8);
  }
  void setupQR() {
    // Set up the Q and R's that we need, apodized.
    const int NkTemp=250;
    std::vector<double> ka(NkTemp),R1(NkTemp),R2(NkTemp);
    std::vector<double> Q1(NkTemp),Q2(NkTemp),Q5(NkTemp),Q8(NkTemp);
    std::vector<double> Qs(NkTemp);
#pragma omp parallel for
    for (int i=1; i<NkTemp; ++i) {
      double kk = exp( kLin[0]+i*(kLin[kLin.size()-1]-kLin[0])/(NkTemp-1) );
      double ap = cos(M_PI/2.*kk/exp(kLin[kLin.size()-1]));
      std::vector<double> Qn=lpt.Qn(kk);
      std::vector<double> Rn=lpt.Rn(kk);
      ka[i]=kk; R1[i]=Rn[0]*ap; R2[i]=Rn[1]*ap;
      Q1[i]=Qn[1]*ap; Q2[i]=Qn[2]*ap; Q5[i]=Qn[5]*ap; Q8[i]=Qn[8]*ap;
      Qs[i]=Qn[14]*ap;
    }
    // and fit splines to them.
    R1spl.init(ka,R1);
    R2spl.init(ka,R2);
    Q1spl.init(ka,Q1);
    Q2spl.init(ka,Q2);
    Q5spl.init(ka,Q5);
    Q8spl.init(ka,Q8);
    Qsspl.init(ka,Qs);
  }
  // This is called something other than xiContributions so that we
  // still have access to the Zeldovich method.
  std::vector<double> dpair(const double rval) {
    // Returns the different contributions to the real-space correlation
    // function for locally biased tracers.
    // This has an abbreviated name so that it doesn't overload the Zeldovich
    // version, allowing comparison and easy switching.
    // This is not tested for very large or small values of r.
    // The integration is over x=q-r, in length and angle with the
    // azimuthal integral being trivial.
    const double xmin=0;
    const double xmax=10*sqrt(sigma2);
    const double rr[3]={0,0,rval};
    const double r2   =rval*rval;
    const int    Nx=500;
    const double dx=(xmax-xmin)/Nx;
    std::vector<double> xi(12);
    for (int ixx=0; ixx<Nx; ++ixx) {
      double xx=xmin+(ixx+0.5)*dx;
      double x2=xx*xx;
      for (int imu=0; imu<gl.N; ++imu) {
        double mu = gl.x[imu];
        // Compute vec{q}=vec{r}+vec{x} with vec{r}=r.zhat,
        // so r_z=r, x_z=x*mu, cos_rq=(r_z+x_z)/q.
        double qlen = sqrt(r2+x2+2*rval*xx*mu);
        double qcos = (rval+xx*mu)/qlen;
        double qsin = sqrt(1-qcos*qcos);
        double qq[3] ={qlen*qsin,0,qlen*qcos};
        double qh[3] ={     qsin,0,     qcos};
        if (qlen>qmin && qlen<qmax) {
          // We keep the Zeldovich piece exponentiated and expand down
          // the 1-loop piece.
          double pref = x2 * zeldovichIntegrand(rr,qq,0) * gl.w[imu];
          // For the bias terms, compute U,xi and Ainv (even though in above).
          std::vector<double> qf  =interpQfuncs(qlen);
          std::vector<double> ef  =interpEfuncs(qlen);
          std::vector<double> Ainv=calcAinv(qq);
          std::vector<double> Aloop(9);
          double Xloop=ef[1]+2*ef[2],Yloop=ef[4]+2*ef[5];
          for (int i=0; i<3; ++i)
            for (int j=i; j<3; ++j) {
              Aloop[3*i+j] = Aloop[3*j+i]  = Yloop*qh[i]*qh[j];
              if (i==j)      Aloop[3*i+i] += Xloop;
            }
          double xiL=qf[3];
          // Construct the auxilliary matrix/vectors g, G of CLPT Eq. (45)
          // and Gamma of Eq. (75).
          double g[3],UL[3],U[3],U20[3],U11[3],G[9];
          for (int i=0; i<3; ++i) {
            g[i]=0;
            for (int j=0; j<3; ++j)
              g[i] += Ainv[3*i+j]*(qq[j]-rr[j]);
            UL[i] = ef[11]*qh[i];
            U[ i] =(ef[11]+ef[12])*qh[i];
            U20[i]=ef[13]*qh[i];
            U11[i]=ef[14]*qh[i];
          }
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              G[3*i+j]=Ainv[3*i+j]-g[i]*g[j];
          double GA=0;
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              GA += Aloop[3*i+j]*G[3*i+j];
          double GW=0;
          double V1=ef[8];
          double V3=ef[9];
          double Tq=ef[10];
          for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
              for (int k=0; k<3; ++k) {
                double Gam,W;
                Gam = Ainv[3*i+j]*g[k]+Ainv[3*k+i]*g[j]+Ainv[3*j+k]*g[i]
                    - g[i]*g[j]*g[k];
                W   = Tq*qh[i]*qh[j]*qh[k];
                if (j==k) W += V1*qh[i];
                if (i==k) W += V1*qh[j];
                if (i==j) W += V3*qh[k];
                GW += Gam*W;
              }
            }
          }
          GW *= 3;	// Account for permutations.
          double trG,Ug,ULg,gq;  trG=Ug=ULg=gq=0;
          for (int i=0; i<3; ++i) {
            gq += g[i]*qh[i];
            ULg+=UL[i]*g[i];
            Ug += U[i]*g[i];
            trG+= G[3*i+i];
          }
          double U20g=0,U11g=0;
          for (int i=0; i<3; ++i) {
            U20g += U20[i]*g[i];
            U11g += U11[i]*g[i];
          }
          double UUG=0,qqG=0;
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j) {
              UUG += G[3*i+j]*UL[i]*UL[j];
              qqG += G[3*i+j]*qh[i]*qh[j];
            }
          double A10G=2*trG*ef[6] + 2*qqG*ef[7];
          double d2xiLin=qf[5];
          // The mode-coupling term, then add the <s^2 Delta Delta> term:
          double shear=ef[15]*gq;
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j) {
              double upsilon= qh[i]*qh[j]*(3*qf[6]*qf[6]+4*qf[6]*qf[7]+
                              2*qf[6]*qf[8]+2*qf[7]*qf[7]+4*qf[7]*qf[8]+
                              qf[8]*qf[8]) + (i==j)*2*qf[7]*qf[7];
              shear += G[3*i+j]*upsilon;
            }
          shear *= 2;
          double V12=qf[9]*gq;
          // Now do the 1, Fp, Fpp, Fp^2, Fp.Fpp, Fpp^2, b_nabla^2,
          // bs, b1.bs2, b2.bs2, bs2^2 terms.
          xi[ 0] +=   pref *(1-GA/2.+GW/6.); // +Aeft*trG*eftNorm);
          xi[ 1] +=  -pref *(2*Ug+A10G);
          xi[ 2] +=  -pref *(UUG+U20g);
          xi[ 3] +=   pref *(xiL-UUG-U11g);
          xi[ 4] +=  -pref *(2*xiL*ULg);
          xi[ 5] +=   pref *xiL*xiL/2;
          xi[ 6] +=  -pref *0.5*trG;
          xi[ 7] +=   pref *d2xiLin;
          xi[ 8] +=  -pref *shear;
          xi[ 9] +=  -pref *2*V12;
          xi[10] +=   pref *qf[10];
          xi[11] +=   pref *qf[11];
        }
      }
    }
    for (int j=0; j<xi.size(); ++j) {
      xi[j] *= dx;      // Convert sum to integral.
      xi[j] *= 2*M_PI;  // The azimuthal integral.
    }
    xi[0] -= 1.0;	// Calculated 1+xi, subtract 1 for xi.
    return(xi);
  }
  // Need to either overload v12, or define something "else".  I
  // choose to define something else, since we may want the ability
  // to switch to v12_Z later for some reason.
  std::vector<double> vpair(const double rval) {
    // Returns the different contributions to the mean infall velocity
    // for locally biased tracers.  Only the line-of-sight component is
    // returned and the result should be multiplied by f and divided by
    // 1+xi(real).
    // This is not tested for very large or small values of r.
    // The integration is over x=q-r, in length and angle with the
    // azimuthal integral being trivial.
    const double xmin=0;
    const double xmax=10*sqrt(sigma2);
    const double rr[3]={0,0,rval};
    const double r2   =rval*rval;
    const int    Nx=500;
    const double dx=(xmax-xmin)/Nx;
    std::vector<double> vv(10);
    for (int ixx=0; ixx<Nx; ++ixx) {
      double xx=xmin+(ixx+0.5)*dx;
      double x2=xx*xx;
      for (int imu=0; imu<gl.N; ++imu) {
        double mu = gl.x[imu];
        // Compute vec{q}=vec{r}+vec{x} with vec{r}=r.zhat,
        // so r_z=r, x_z=x*mu, cos_rq=(r_z+x_z)/q.
        double qlen = sqrt(r2+x2+2*rval*xx*mu);
        double qcos = (rval+xx*mu)/qlen;
        double qsin = sqrt(1-qcos*qcos);
        double qq[3] ={qlen*qsin,0,qlen*qcos};
        double qh[3] ={     qsin,0,     qcos};
        if (qlen>qmin && qlen<qmax) {
          // We keep the Zeldovich piece exponentiated and expand down
          // the 1-loop piece.
          double pref = x2 * zeldovichIntegrand(rr,qq,0) * gl.w[imu];
          // For the bias terms, compute U,xi and Ainv (even though in above).
          std::vector<double> qf  =interpQfuncs(qlen);
          std::vector<double> ef  =interpEfuncs(qlen);
          std::vector<double> Ainv=calcAinv(qq);
          std::vector<double> Alin(9);
          std::vector<double> Adot(9);
          std::vector<double> Aloop(9);
          double Xdot=ef[0]+2*ef[1]+4*ef[2],Ydot=ef[3]+2*ef[4]+4*ef[5];
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j) {
              Aloop[3*i+j] = (ef[4]+2*ef[5])*qh[i]*qh[j]+(ef[1]+2*ef[2])*(i==j);
              Adot[ 3*i+j] = Ydot*qh[i]*qh[j]+Xdot*(i==j);
              Alin[ 3*i+j] = ef[3]*qh[i]*qh[j]+ef[0]*(i==j);
            }
          double xiL=qf[3];
          // Construct the auxilliary matrix/vectors g, G of CLPT Eq. (45)
          // and Gamma of Eq. (75).
          double g[3],UL[3],U[3],Udot[3],U20[3],U11[3],G[9];
          for (int i=0; i<3; ++i) {
            g[i]=0;
            for (int j=0; j<3; ++j)
              g[i] += Ainv[3*i+j]*(qq[j]-rr[j]);
            UL[i]    = ef[11]*qh[i];
            U[ i]    =(ef[11]+  ef[12])*qh[i];
            Udot[i]  =(ef[11]+3*ef[12])*qh[i];
            U20[i]   = 2*ef[13]*qh[i];
            U11[i]   = 2*ef[14]*qh[i];
          }
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              G[3*i+j]=Ainv[3*i+j]-g[i]*g[j];
          double GA=0;
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              GA += Aloop[3*i+j]*G[3*i+j];
          double trG,Ug,gq,qG,gA,UGA,qGq;  trG=Ug=gq=qG=gA=UGA=qGq=0;
          double ULg=0,gAL=0;
          for (int i=0; i<3; ++i) {
            gq += g[i]*qh[i];
            gA += g[i]*Adot[3*2+i];
            gAL+= g[i]*Alin[3*2+i];
            Ug += U[i]*g[i];
            ULg+=UL[i]*g[i];
            qG += qh[i]*G[3*2+i];
            trG+= G[3*i+i];
            for (int j=0; j<3; ++j) {
              UGA += UL[i]*G[3*i+j]*Alin[3*2+j];
              qGq += qh[i]*G[3*i+j]*qh[j];
            }
          }
          double gA10=3*(ef[6]*g[2]+ef[7]*gq*qh[2]);
          double V1=ef[8];
          double V3=ef[9];
          double Tq=ef[10];
          double GW=0;
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j) {
              double Wdot_ijn=(3*V1+V3)*(qh[i]*(j==2)+qh[j]*(i==2))+
                              2*(V1+V3)*qh[2]*(i==j)+
                              4*Tq*qh[i]*qh[j]*qh[2];
              GW += G[3*i+j]*Wdot_ijn;
            }
          // The mode-coupling term, then add the <s^2 Delta Delta> term:
          double shear=2*ef[15]*g[2];
          for (int i=0; i<3; ++i) {
            double upsilon= qh[i]*qh[2]*(3*qf[6]*qf[6]+4*qf[6]*qf[7]+
                            2*qf[6]*qf[8]+2*qf[7]*qf[7]+4*qf[7]*qf[8]+
                            qf[8]*qf[8]) + (i==2)*2*qf[7]*qf[7];
            shear -= g[i]*upsilon;
          }
          shear *= 2;
          double V12=1*qf[9]*gq;
          // Now do the 1, Fp, Fpp, Fp^2, Fp.Fpp, Fpp^2, grad_xi, g_los,
          // bs2 terms.
          vv[0] +=   -pref *(gA+0.5*GW);
          vv[1] +=  2*pref *(Udot[2]-UGA-gA10);
          vv[2] +=    pref *(U20[2]-2*ULg*UL[2]);
          vv[3] +=   -pref *(ULg*UL[2]+ULg*UL[2]+xiL*gAL-U11[2]);
          vv[4] +=  2*pref *xiL*UL[2];
          vv[5] +=  0;
          vv[6] +=    pref *qf[4];
          vv[7] +=    pref *g[2];
          vv[8] +=    pref *shear;
          vv[9] +=    pref *V12;
        }
      }
    }
    for (int j=0; j<vv.size(); ++j) {
      vv[j] *= dx;      // Convert sum to integral.
      vv[j] *= 2*M_PI;  // The azimuthal integral.
    }
    return(vv);
  }
  // Need to either overload s12, or define something "else".  I
  // choose to define something else, since we may want the ability
  // to switch to s12_Z later for some reason.
  std::vector<double> spair(const double rval) {
    // Returns the different contributions to the velocity dispersion
    // for locally biased tracers.  Both sigma_perp and sigma_par are
    // returned and the result should be multiplied by f^2 and divided by
    // 1+xi(real).  Beware the ordering here!!
    // This is not tested for very large or small values of r.
    // The integration is over x=q-r, in length and angle with the
    // azimuthal integral being trivial.
    const double xmin=0;
    const double xmax=10*sqrt(sigma2);
    const double rr[3]={0,0,rval};
    const double r2   =rval*rval;
    const int    Nx=500;
    const double dx=(xmax-xmin)/Nx;
    std::vector<double> ss(16);
    for (int ixx=0; ixx<Nx; ++ixx) {
      double xx=xmin+(ixx+0.5)*dx;
      double x2=xx*xx;
      for (int imu=0; imu<gl.N; ++imu) {
        double mu = gl.x[imu];
        // Compute vec{q}=vec{r}+vec{x} with vec{r}=r.zhat,
        // so r_z=r, x_z=x*mu, cos_rq=(r_z+x_z)/q.
        double qlen = sqrt(r2+x2+2*rval*xx*mu);
        double qcos = (rval+xx*mu)/qlen;
        double qsin = sqrt(1-qcos*qcos);
        double qq[3] ={qlen*qsin,0,qlen*qcos};
        double qh[3] ={     qsin,0,     qcos};
        if (qlen>qmin && qlen<qmax) {
          // For the unbiased tracers we only need this--all other terms
          // are multiplied by this anyway.
          double pref = x2 * zeldovichIntegrand(rr,qq,0) * gl.w[imu];
          // For the bias terms compute U, xi and Ainv (even though in above).
          std::vector<double> qf  =interpQfuncs(qlen);
          std::vector<double> ef  =interpEfuncs(qlen);
          std::vector<double> Ainv=calcAinv(qq);
          std::vector<double> Alin(9);
          std::vector<double> Adot(9);
	  std::vector<double> Addot(9);
          std::vector<double> Aloop(9);
          double Xdot =ef[0]+2*ef[1]+4*ef[2],Ydot =ef[3]+2*ef[4]+4*ef[5];
          double Xddot=ef[0]+4*ef[1]+6*ef[2],Yddot=ef[3]+4*ef[4]+6*ef[5];
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j) {
              Aloop[3*i+j] = (ef[4]+2*ef[5])*qh[i]*qh[j]+(ef[1]+2*ef[2])*(i==j);
              Alin[ 3*i+j] = ef[3]*qh[i]*qh[j]+ef[0]*(i==j);
              Adot[ 3*i+j] =  Ydot*qh[i]*qh[j]+ Xdot*(i==j);
              Addot[3*i+j] = Yddot*qh[i]*qh[j]+Xddot*(i==j);
            }
          double xiL=qf[3];
          // Construct the auxilliary matrix/vectors g, G of CLPT Eq. (45)
          double g[3],UL[3],U[3],Udot[3],G[9],W[27],Wddot[27];
          for (int i=0; i<3; ++i) {
            g[i]=0;
            for (int j=0; j<3; ++j)
              g[i] += Ainv[3*i+j]*(qq[j]-rr[j]);
            UL[i]   = ef[11]*qh[i];
            U[ i]   =(ef[11]+  ef[12])*qh[i];
            Udot[i] =(ef[11]+3*ef[12])*qh[i];
          }
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              G[3*i+j]=Ainv[3*i+j]-g[i]*g[j];
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              for (int k=0; k<3; ++k)
                W[9*i+3*j+k] = ef[8]*qh[i]*(j==k) + ef[8]*qh[j]*(i==k)
                             + ef[9]*qh[k]*(i==j) +ef[10]*qh[i]*qh[j]*qh[k];
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              for (int k=0; k<3; ++k)
                Wddot[9*i+3*j+k] = 2*W[9*i+3*j+k]+2*W[9*i+3*k+j]+W[9*k+3*j+i];
          // We also need \ddot{A}^{10}:
          double A10[9];
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
                A10[3*i+j] = (4*ef[7])*qh[i]*qh[j] + (4*ef[6])*(i==j);
          // Work out these dot products for n=m=2, i.e. along z, for sigma_par.
          double Ug,ULg,U2,gA,gAL,gAU,AGA;  Ug=ULg=U2=gA=gAL=gAU=AGA=0;
          for (int i=0; i<3; ++i) {
            gA += g[i]*Adot[3*2+i];
            gAL+= g[i]*Alin[3*2+i];
            Ug += U[i]*g[i];
            ULg+=UL[i]*g[i];
            U2 +=UL[i]*UL[i];
            for (int j=0; j<3; ++j) {
              gAU += g[i]*Alin[3*i+j]*UL[j];
              AGA += Alin[3*2+i]*G[3*i+j]*Alin[3*2+j];
            }
          }
          double Wg=0;
          for (int i=0; i<3; ++i) Wg += Wddot[9*i+3*2+2]*g[i];
          // Now the shear term.
          double upsilon= qh[2]*qh[2]*(3*qf[6]*qf[6]+4*qf[6]*qf[7]+
                          2*qf[6]*qf[8]+2*qf[7]*qf[7]+4*qf[7]*qf[8]+
                          qf[8]*qf[8]) + 2*qf[7]*qf[7];
          double shear = 2*upsilon;
          // Now do the 1, Fp, Fpp, Fp^2, Fp.Fpp, Fpp^2 terms for \sigma_par^2.
          ss[ 0] +=    pref *(Addot[3*2+2]-AGA-Wg);
          ss[ 1] += -2*pref *(2*gAL*UL[2]+ULg*Alin[3*2+2]-A10[3*2+2]);
          ss[ 2] +=  2*pref *   UL[2]*UL[2];
          ss[ 3] +=    pref *(2*UL[2]*UL[2]+xiL*Alin[3*2+2]);
          ss[ 4]  =  0;
          ss[ 5]  =  0;
          ss[ 6]  =    pref *2*shear;
          ss[ 7]  =    pref *1*qf[3];
          // Next work out the trace components, i.e. summed over n=m.
          Wg=0;
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              Wg += Wddot[9*i+3*j+j]*g[i];
          AGA=0;
          for (int m=0; m<3; ++m) {
            for (int i=0; i<3; ++i) {
              for (int j=0; j<3; ++j) {
                AGA += Alin[3*m+i]*G[3*i+j]*Alin[3*j+m];
              }
            }
          }
          double trA=0,trAL=0,trA10=0;
          for (int m=0; m<3; ++m) {
            trA  += Addot[3*m+m];
            trAL +=  Alin[3*m+m];
            trA10+=   A10[3*m+m];
          }
          // Now the shear term.
          upsilon=0;
          for (int m=0; m<3; ++m)
            upsilon += qh[m]*qh[m]*(3*qf[6]*qf[6]+4*qf[6]*qf[7]+
                       2*qf[6]*qf[8]+2*qf[7]*qf[7]+4*qf[7]*qf[8]+
                       qf[8]*qf[8]) + 2*qf[7]*qf[7];
          shear = 2*upsilon;
          ss[ 8] +=    pref *(trA-AGA-Wg);
          ss[ 9] += -2*pref *(2*gAU+ULg*trAL-trA10);
          ss[10] +=  2*pref * U2;
          ss[11] +=    pref *(2*U2+xiL*trAL);
          ss[12]  =  0;
          ss[13]  =  0;
          ss[14]  =    pref *2*shear;
          ss[15]  =    pref *3*qf[3];
        }
      }
    }
    // Now sigma_perp is related to the trace by s_perp^2=(1/2)[Tr-sig_par^2]
    for (int j=8; j<ss.size(); ++j)
      ss[j] = 0.5*(ss[j] - ss[j-8]);
    for (int j=0; j<ss.size(); ++j) {
      ss[j] *= dx;	// Convert sum to integral.
      ss[j] *= 2*M_PI;	// The azimuthal integral.
    }
    // Note we return parallel then perpendicular/transverse!
    return(ss);
  }
public:
  LSM() {}	// Do nothing.
  LSM(const char fname[], const double f,
      const double b1, const double b2, const double bs,
      const double Aeft, const double Aeft1, const double Aeft2) {
    init(fname,f,b1,b2,bs,Aeft,Aeft1,Aeft2);
  }
  void init(const char fname[], const double f,
            const double b1, const double b2, const double bs,
            const double Aeft, const double Aeft1, const double Aeft2) {
    // Set up the Zel'dovich class.
    Zeldovich::init(fname);
    // Initialize the LPT class with our newly populated kLin/pLin.
    lpt.init(kLin,pLin);
    // We can now set up the "extra" functions we need.  Check
    // to see whether we have this pretabulated.
    std::stringstream ss;
    ss<<fname<<".lesmSave";
    std::ifstream fs(ss.str().c_str());
    if (!fs) {
      // Set up the Q and R's that we need.
      setupQR();
      // and tabulate/calculate the "extra" functions.
      tabulateEfuncs();
      writeSaveFile(ss.str().c_str());
    }
    else {
      fs.close();
      readSaveFile(ss.str().c_str());
    }
    // Now tabulate the functions we need for the streaming model.
    // The point at zero lag is known analytically.
    std::vector<double>	rrvec,xivec,vvvec,stvec,spvec;
    double rr,xi,zi,vv,st,sp;
    // Step up to the maximum distance.
    rr=xi=vv=st=sp=0;
    const double dr=2,rmax=250;
    do {
      try {
        rrvec.push_back(rr);
      } catch(std::exception& e) {myexception(e);}
      rr += dr;
    } while(rr<rmax);
    try {
      xivec.resize( rrvec.size() );
      vvvec.resize( rrvec.size() );
      stvec.resize( rrvec.size() );
      spvec.resize( rrvec.size() );
    } catch(std::exception& e) {myexception(e);}
#pragma omp parallel for
    for (int i=1; i<rrvec.size(); ++i) {
      //std::vector<double> zC = xiContributions(rrvec[i],Aeft);
      //std::vector<double> vC = v12(rrvec[i]);
      std::vector<double> xC = dpair(rrvec[i]);
      std::vector<double> vC = vpair(rrvec[i]);
      std::vector<double> sC = spair(rrvec[i]);
      xi = xC[0]+b1*xC[1]+b2*xC[2]+b1*b1*xC[3]+b1*b2*xC[ 4]+b2*b2*xC[ 5]+
               Aeft*xC[6]+ 0*xC[7]+   bs*xC[8]+b1*bs*xC[ 9]+b2*bs*xC[10]+
              bs*bs*xC[11];
      vv = vC[0]+b1*vC[1]+b2*vC[2]+b1*b1*vC[3]+b1*b2*vC[ 4]+b2*b2*vC[ 5]+
              Aeft1*vC[6]+Aeft2*vC[7]+bs*vC[8]+b1*bs*vC[ 9];
      sp = sC[0]+b1*sC[1]+b2*sC[2]+b1*b1*sC[ 3]+b1*b2*sC[ 4]+b2*b2*sC[ 5]+
           bs*sC[6]+0*sC[7];
      st = sC[8]+b1*sC[9]+b2*sC[10]+b1*b1*sC[11]+b1*b2*sC[12]+b2*b2*sC[13]+
           bs*sC[14]+0*sC[15];
      vv*=   f/(1+xi);
      sp*= f*f/(1+xi);
      st*= f*f/(1+xi);
      xivec[i] = xi*rrvec[i]*rrvec[i];  // Actually stores r^2.xi
      vvvec[i] = vv;
      spvec[i] = sp;
      stvec[i] = st;
    }
    // and fit splines to them for later use.
    xispl.init(rrvec,xivec);
    vvspl.init(rrvec,vvvec);
    stspl.init(rrvec,stvec);
    spspl.init(rrvec,spvec);
  }
  double xiRZ(const double R, const double Z, const double s2fog) {
    // The 2D correlation function for the streaming model.
    // Does the integral over the "true" line-of-sight separation
    // using Simpson's rule.
    const double R2=R*R;
    const double ymax=50;
    const int    Ny=500;
    const double hh=2*ymax/Ny;
    int errcnt=0;
    double xi=0;
    // Careful throwing exceptions from threads...
#pragma omp parallel for reduction(+:xi,errcnt)
    for (int i=1; i<Ny; ++i) {
      double yy  = -ymax + i*hh;        // Actually Z-y
      double rr  = sqrt(R2+(Z-yy)*(Z-yy));
      if (errcnt>0 || rr>=xispl.xmax() || rr<=xispl.xmin()) {
        errcnt=1;
      }
      else {
        double mu  = (Z-yy)/rr;
        double xip1= 1.0 + xispl(rr)/rr/rr;
        double vr  = mu*vvspl(rr);
        double expt= yy-vr;
        double s2  = mu*mu*spspl(rr)+(1-mu*mu)*stspl(rr)-vr*vr + s2fog;
        int    wt  = 2+2*(i%2);
        if (s2>0)
          xi += xip1*exp(-0.5*expt*expt/s2)/sqrt(s2) * wt;
      }
    }
    if (errcnt>0) {myexit(1);}
    xi *= hh/3.0 / sqrt(2*M_PI);
    xi -= 1.0;
    return(xi);
  }
  std::vector<double> xiEll(const double ss, const double s2fog,
                            const double Apar, const double Aperp) {
    // The multipoles of the correlation function for the streaming model.
    // Integrates the 2D correlation function using Gauss-Legendre integration.
    std::vector<double> xiell(2);
    const int Nmu=16;
    GaussLegendre gg = GaussLegendre(2*Nmu);	// Must be even.
    // For even lengths, can sum over half of the points.
    for (int i=0; i<Nmu; ++i) {
      double ximu = xiRZ(ss*sqrt(1-gg.x[i]*gg.x[i])*Aperp,
                         ss*gg.x[i]*Apar,s2fog);
      double p0=1.0;
      double p2=0.5*(3*gg.x[i]*gg.x[i]-1);
      xiell[0] += ximu*gg.w[i] * p0 * 1;
      xiell[1] += ximu*gg.w[i] * p2 * 5;
    }
    return(xiell);
  }
  void printzFuncs(const char fbase[]) {
    // Print the "extra" functions.
    std::ostringstream ss;
    ss << fbase << ".zFuncs";
    std::ofstream fs(ss.str().c_str());
    if (!fs) {std::cerr<<"Unable to open file."<<std::endl;myexit(1);}
    fs<<"# q-dependent functions computed for Zeldovich."<<std::endl;
    fs<<"# Order is"<<std::endl;
    fs<<"#  1) q [Mpc/h]"<<std::endl
      <<"#  2) eta_per" <<std::endl
      <<"#  3) eta_par" <<std::endl
      <<"#  4) U^{(1)}" <<std::endl
      <<"#  5) xi_L"    <<std::endl
      <<"#  6) xi_L\'"  <<std::endl
      <<"#  7) d^2 xi_L"<<std::endl
      <<"#  8) J_2"     <<std::endl
      <<"#  9) J_3"     <<std::endl
      <<"# 10) J_4"     <<std::endl
      <<"# 11) V_i^{12}"<<std::endl
      <<"# 12) chi^{12}"<<std::endl
      <<"# 13) zeta"    <<std::endl;
    for (int i=0; i<120; ++i) {
      double q = (i+1.0);
      std::vector<double> qf = interpQfuncs(q);
      fs<<std::scientific<<std::setw(15)<<std::setprecision(5)<<q;
      for (int j=0; j<qf.size(); ++j)
        fs<<std::scientific<<std::setw(15)<<std::setprecision(5)<<qf[j];
      fs<<std::endl;
    }
    fs.close();
  }
  void printqFuncs(const char fbase[]) {
    // Print the "extra" functions.
    std::ostringstream ss;
    ss << fbase << ".qFuncs";
    std::ofstream fs(ss.str().c_str());
    if (!fs) {std::cerr<<"Unable to open file."<<std::endl;myexit(1);}
    fs<<"# q-dependent functions used by CLPT."<<std::endl;
    fs<<"# Order is"<<std::endl;
    fs<<"#  1) q [Mpc/h]"<<std::endl
      <<"#  2) X^{(11)}"<<std::endl
      <<"#  3) X^{(22)}"<<std::endl
      <<"#  4) X^{(13)}"<<std::endl
      <<"#  5) Y^{(11)}"<<std::endl
      <<"#  6) Y^{(22)}"<<std::endl
      <<"#  7) Y^{(13)}"<<std::endl
      <<"#  8) X^{(12)}_{10}"<<std::endl
      <<"#  9) Y^{(12)}_{10}"<<std::endl
      <<"# 10) V^{(112)}_{1}"<<std::endl
      <<"# 11) V^{(112)}_{3}"<<std::endl
      <<"# 12) T^{(112)}"<<std::endl
      <<"# 13) U^{(1)}"<<std::endl
      <<"# 14) U^{(3)}"<<std::endl
      <<"# 15) U^{(2)}_{20}"<<std::endl
      <<"# 16) U^{(2)}_{11}"<<std::endl
      <<"# 17) V^{10}"<<std::endl;
    for (int i=0; i<60; ++i) {
      double q = 2*(i+0.5);
      std::vector<double> ef = interpEfuncs(q);
      fs<<std::scientific<<std::setw(15)<<std::setprecision(5)<<q;
      for (int j=0; j<ef.size(); ++j)
        fs<<std::scientific<<std::setw(15)<<std::setprecision(5)<<ef[j];
      fs<<std::endl;
    }
    fs.close();
  }
  void printXiStuff(const char fbase[]) {
    // Print the contributions to Xi
    std::ostringstream ss;
    ss << fbase << ".xiStuff";
    std::ofstream fs(ss.str().c_str());
    if (!fs) {std::cerr<<"Unable to open file."<<std::endl;myexit(1);}
    fs << "# Contributions to xi_real."<<std::endl
       << "# Order is r [Mpc/h], xi_L,"
       << " 1, b1, b2, b1^2, b1.b2, b2^2, Aeft, d2xiLin, bs, b1.bs, b2.bs, bs^2"
       << std::endl;
    for (int ir=0; ir<65; ++ir) {
      double rr = 10.0 + 2*(ir+0.5);
      std::vector<double> qf = interpQfuncs(rr);
      std::vector<double> xC = dpair(rr);
      fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<rr
        <<std::scientific<<std::setw(12)<<std::setprecision(4)<<qf[3];
      for (int j=0; j<xC.size(); ++j)
        fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<xC[j];
      fs<<std::endl;
    }
    fs.close();
  }
  void printVpStuff(const char fbase[]) {
    // Print the contributions to Vpair
    std::ostringstream ss;
    ss << fbase << ".vpStuff";
    std::ofstream fs(ss.str().c_str());
    if (!fs) {std::cerr<<"Unable to open file."<<std::endl;myexit(1);}
    fs << "# Contributions to (1/f)(1+xi)v_{12}."<<std::endl
       << "# Order is r [Mpc/h] "
       << "v_L, 1, b1, b2, b1^2, b1.b2, b2^2, gradXi, g_los, bs2, bs2.b1"
       << std::endl;
    for (int ir=0; ir<65; ++ir) {
      double rr = 10.0 + 2*(ir+0.5);
      std::vector<double> qf = interpQfuncs(rr);
      std::vector<double> vC = vpair(rr);
      fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<rr
        <<std::scientific<<std::setw(12)<<std::setprecision(4)<<2*qf[2];
      for (int j=0; j<vC.size(); ++j)
        fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<vC[j];
      fs<<std::endl;
    }
    fs.close();
  }
  void printS2Stuff(const char fbase[]) {
    // Print the contributions to Sigma^2
    std::ostringstream ss;
    ss << fbase << ".s2Stuff";
    std::ofstream fs(ss.str().c_str());
    if (!fs) {std::cerr<<"Unable to open file."<<std::endl;myexit(1);}
    fs << "# Contributions to (1/f^2)(1+xi)sigma^2."<<std::endl
       << "# Order is r [Mpc/h]"<<std::endl
       << "# 1, b1, b2, b1^2, b1.b2, b2^2, bs2, beta"<<std::endl
       << "# first for sig2_par then for sig2_perp."<<std::endl;
    for (int ir=0; ir<65; ++ir) {
      double rr = 10.0 + 2*(ir+0.5);
      std::vector<double> sC = spair(rr);
      fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<rr;
      for (int j=0; j<sC.size(); ++j)
        fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<sC[j];
      fs<<std::endl;
    }
    fs.close();
  }
};




int	main(int argc, char **argv)
{
  if (argc != 9) {
    std::cout<<"Usage: lesm "
             <<"<Pk-file> <ff> <b1> <b2> <bs2> <Aeft> <Aeftv> <s2FoG>"
             <<std::endl;
    myexit(1);
  }
  const double ff   = atof(argv[2]);
  const double b1   = atof(argv[3]);
  const double b2   = atof(argv[4]);
  const double bs2  = atof(argv[5]);
  const double Apar = 1;
  const double Aperp= 1;
  const double Aeft = atof(argv[6]);
  const double Aeft1= atof(argv[7]);
  const double Aeft2= 0;
  const double s2FoG= atof(argv[8]);
  try {
    LSM lsm(argv[1],ff,b1,b2,bs2,Aeft,Aeft1,Aeft2);

#ifdef	PRINTSTUFF
    lsm.printzFuncs( argv[1]);
    lsm.printqFuncs( argv[1]);
    lsm.printXiStuff(argv[1]);
    lsm.printVpStuff(argv[1]);
    lsm.printS2Stuff(argv[1]);
    return(0);
#endif

    const int Nr=61;
    const double rmin=5.0,rmax=135;
    std::vector<double> xiell;
    for (int i=0; i<Nr; ++i) {
      double rr = rmin + i*(rmax-rmin)/(Nr-1);
      xiell=lsm.xiEll(rr,s2FoG,Apar,Aperp);
      std::cout<<std::fixed<<std::setw(10)<<std::setprecision(2)<<rr
               <<std::fixed<<std::setw(15)<<std::setprecision(5)<<xiell[0]*rr*rr
               <<std::fixed<<std::setw(15)<<std::setprecision(5)<<xiell[1]*rr*rr
               <<std::endl;
    }
  } catch(std::exception& e) {myexception(e);}
  return(0);
}
