#include <cmath>

#include "lpt.hpp"

double
LPT::linInterp(const double kk) const {
    // Returns a linear interpolant to P(k).
    double lnk = log(kk);
    if (lnk<=kLin[0] || lnk>=kLin[kLin.size()-1]) {return(0);}
    int jj = (int)((lnk-kLin[0])*dkinv);
    if (jj>=pLin.size()-2) jj=pLin.size()-2;
    double pk  = exp(pLin[jj]+(lnk-kLin[jj])*(pLin[jj+1]-pLin[jj])/(kLin[jj+1]-kLin[jj]));
    return(pk);
}

std::vector<double>
LPT::Qdefs(const double r, const double x) const {
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

std::vector<double>
LPT::Rdefs(const double r) const {
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

std::vector<double>
LPT::Qnx(const double kval, const double rr) const {
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

std::vector<double>
LPT::Qn(const double kk) const {
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

std::vector<double>
LPT::Rn(const double kk) const {
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

void
LPT::init(std::vector<double> inkLin, std::vector<double> inpLin) {
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
