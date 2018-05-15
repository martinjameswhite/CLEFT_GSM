#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "utils.hpp"
#include "spline.hpp"
#include "zeldovich.hpp"
#include "four1.hpp"

std::vector<double>
Zeldovich::sphBess(const double x) {
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

void
Zeldovich::readPowerSpectrum(const char fname[]) {
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

double
Zeldovich::calcSigma2() {
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

double
Zeldovich::calcEFTnorm() {
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

std::vector<double>
Zeldovich::calcQfuncs(const double q) {
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

double
Zeldovich::calc_nabla1(const double q) {
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

double
Zeldovich::calc_nabla2(const double q) {
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

std::vector<double>
Zeldovich::calc_Jn(const double q) {
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

void
Zeldovich::tabulateQfuncs() {
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

std::vector<double>
Zeldovich::interpQfuncs(const double q) {
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

std::vector<double>
Zeldovich::calcAmat(const double q[]) {
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

std::vector<double>
Zeldovich::calcAinv(const double q[]) {
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

double
Zeldovich::zeldovichIntegrand(const double r[], const double q[], const double f){
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

void
Zeldovich::init(const char fname[]) {
    // Initialize the G-L integration points and weights.
    gl.set(128);
    // Load the linear power spectrum from a file, and expand it out
    // if necessary.
    readPowerSpectrum(fname);
    // Finally pre-tabulate the q-dependent functions.
    tabulateQfuncs();
    // and the normalization term for EFT.
    eftNorm = calcEFTnorm();
}

void
Zeldovich::print_eta() {
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

double
Zeldovich::xiL(const double r) {
    // The real-space, linear correlation function at r.
    // This is not tested for very large or small values of r.
    std::vector<double> qf=interpQfuncs(r);
    return(qf[3]);
}

double
Zeldovich::xiZ(const double rval) {
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

std::vector<double>
Zeldovich::xiContributions(const double rval, const double Aeft) {
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

std::vector<double>
Zeldovich::v12(const double rval) {
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

std::vector<double>
Zeldovich::s12(const double rval) {
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
