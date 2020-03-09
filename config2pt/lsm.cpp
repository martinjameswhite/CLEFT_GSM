#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "utils.hpp"
#include "lsm.hpp"

std::vector<double>
LSM::calcEfuncs(const double q) {
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
      sum12+= k2*wt*(Qs*j1)*(-1./7.)*ap;		// Shear term
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

void
LSM::tabulateEfuncs() {
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

void
LSM::writeSaveFile(const char fname[]) {
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

void
LSM::readSaveFile(const char fname[]) {
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

std::vector<double>
LSM::interpEfuncs(const double q) {
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

void
LSM::setupQR() {
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


// Returns the different contributions to the real-space correlation function (component 0)
// mean infall velocity (component 1), and the velocity dispersion (component 2)
// for locally biased tracers.
// For mean infall velocity (component 1) only the line-of-sight component is returned
// and the result should be multiplied by f and divided by 1+xi(real).
// For velocity dispersion (component 2) both sigma_perp and sigma_par are returned
// and the result should be multiplied by f^2 and divided by 1+xi(real). NOTE we return
// parallel then perpendicular/transverse.
// This is not tested for very large or small values of r.
// The integration is over x=q-r, in length and angle with the
// azimuthal integral being trivial.
std::vector<std::vector<double>>
LSM::dvsPair(const double rval)
{
    const double xmin=0;
    const double xmax=10*sqrt(sigma2);
    const double rr[3]={0,0,rval};
    const double r2   =rval*rval;
    const int    Nx=500;
    const double dx=(xmax-xmin)/Nx;

    std::vector<double> xi(12);
    std::vector<double> vv(10);
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
            double qq[3] = {qlen*qsin,0,qlen*qcos};
            double qh[3] = {     qsin,0,     qcos};
            if (qlen>qmin && qlen<qmax) {
                // We keep the Zeldovich piece exponentiated and expand down
                // the 1-loop piece.
                double pref = x2 * zeldovichIntegrand(rr,qq,0) * gl.w[imu];
                // For the bias terms, compute U, xi and Ainv (even though in above).
                std::vector<double> qf  =interpQfuncs(qlen);
                std::vector<double> ef  =interpEfuncs(qlen);
                std::vector<double> Ainv=calcAinv(qq);
                std::vector<double> Aloop(9);
                std::vector<double> Alin(9);
                std::vector<double> Adot(9);
                std::vector<double> Addot(9);

                double Xdot=ef[0]+2*ef[1]+4*ef[2],Ydot=ef[3]+2*ef[4]+4*ef[5];
                double Xddot=ef[0]+4*ef[1]+6*ef[2],Yddot=ef[3]+4*ef[4]+6*ef[5];
                for (int i=0; i<3; ++i) {
                    for (int j=0; j<3; ++j) {
                        Aloop[3*i+j] = (ef[4]+2*ef[5])*qh[i]*qh[j]+(ef[1]+2*ef[2])*(i==j);
                        Adot[ 3*i+j] = Ydot*qh[i]*qh[j]+Xdot*(i==j);
                        Alin[ 3*i+j] = ef[3]*qh[i]*qh[j]+ef[0]*(i==j);
                        Addot[3*i+j] = Yddot*qh[i]*qh[j]+Xddot*(i==j);
                    }
                }
                double xiL=qf[3];
                // Construct the auxilliary matrix/vectors g, G of CLPT Eq. (45)
                // and Gamma of Eq. (75).
                double g[3],UL[3],U[3],Udot[3],U20[3],U11[3],G[9],W[27],Wddot[27];
                for (int i=0; i<3; ++i) {
                    g[i]=0;
                    for (int j=0; j<3; ++j)
                        g[i] += Ainv[3*i+j]*(qq[j]-rr[j]);
                    UL[i] = ef[11]*qh[i];
                    U[ i] =(ef[11]+ef[12])*qh[i];
                    Udot[i]  =(ef[11]+3*ef[12])*qh[i];
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
                double trG = 0, Ug = 0, ULg = 0, U2 = 0, gq = 0, qG = 0, gA = 0, UGA = 0, qGq = 0, gAL = 0, gAU = 0, AGA = 0;
                for (int i=0; i<3; ++i) {
                    gq += g[i]*qh[i];
                    gA += g[i]*Adot[3*2+i];
                    gAL+= g[i]*Alin[3*2+i];
                    Ug += U[i]*g[i];
                    ULg+=UL[i]*g[i];
                    U2 +=UL[i]*UL[i];
                    qG += qh[i]*G[3*2+i];
                    trG+= G[3*i+i];
                    for (int j=0; j<3; ++j) {
                        UGA += UL[i]*G[3*i+j]*Alin[3*2+j];
                        qGq += qh[i]*G[3*i+j]*qh[j];
                        gAU += g[i]*Alin[3*i+j]*UL[j];
                        AGA += Alin[3*2+i]*G[3*i+j]*Alin[3*2+j];
                    }
                }

                double GWv=0;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j) {
                        double Wdot_ijn=(3*V1+V3)*(qh[i]*(j==2)+qh[j]*(i==2))+
                                  2*(V1+V3)*qh[2]*(i==j)+
                                  4*Tq*qh[i]*qh[j]*qh[2];
                        GWv += G[3*i+j]*Wdot_ijn;
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

                double gA10=3*(ef[6]*g[2]+ef[7]*gq*qh[2]);

                // The mode-coupling term, then add the <s^2 Delta Delta> term:
                double shear=ef[15]*gq;
                for (int i=0; i<3; ++i)
                    for (int j=0; j<3; ++j) {
                        double upsilon2= qh[i]*qh[j]*(3*qf[6]*qf[6]+
                                  4*qf[6]*qf[7]+2*qf[6]*qf[8]+2*qf[7]*qf[7]+
                                  4*qf[7]*qf[8]+  qf[8]*qf[8]) +
                                  (i==j)*2*qf[7]*qf[7];
                        shear += G[3*i+j]*upsilon2;
                    }
                shear *= 2;
                // The mode-coupling term, then add the <s^2 Delta Delta> term:
                double shear_v=4*ef[15]*qh[2];
                for (int i=0; i<3; ++i) {
                    double upsilon2= qh[i]*qh[2]*(3*qf[6]*qf[6]+4*qf[6]*qf[7]+
                                    2*qf[6]*qf[8]+2*qf[7]*qf[7]+4*qf[7]*qf[8]+
                                      qf[8]*qf[8]) + (i==2)*2*qf[7]*qf[7];
                    shear_v -= 4*g[i]*upsilon2;
                }
                double V12 =qf[9]*gq;
                double V12z=qf[9]*qh[2];

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

                // Now do the 1, Fp, Fpp, Fp^2, Fp.Fpp, Fpp^2, grad_xi, g_los,
                // bs2 terms.
                vv[0] +=   -pref *(gA+0.5*GWv);
                vv[1] +=  2*pref *(Udot[2]-UGA-gA10);
                vv[2] +=    pref *(2*U20[2]-2*ULg*UL[2]);
                vv[3] +=   -pref *(ULg*UL[2]+ULg*UL[2]+xiL*gAL-2*U11[2]);
                vv[4] +=  2*pref *xiL*UL[2];
                vv[5] +=  0;
                vv[6] +=    pref *qf[4];
                vv[7] +=    pref *g[2];
                vv[8] +=    pref *shear_v;
                vv[9] +=  2*pref *V12z;

                double Wg=0;
                for (int i=0; i<3; ++i) Wg += Wddot[9*i+3*2+2]*g[i];
                // Now the shear term.
                double upsilon= qh[2]*qh[2]*(3*qf[6]*qf[6]+4*qf[6]*qf[7]+
                          2*qf[6]*qf[8]+2*qf[7]*qf[7]+4*qf[7]*qf[8]+
                          qf[8]*qf[8]) + 2*qf[7]*qf[7];
                double shear_s = 2*upsilon;
                // Now do the 1, Fp, Fpp, Fp^2, Fp.Fpp, Fpp^2 terms for \sigma_par^2.
                ss[ 0] +=    pref *(Addot[3*2+2]-AGA-Wg);
                ss[ 1] += -2*pref *(2*gAL*UL[2]+ULg*Alin[3*2+2]-A10[3*2+2]);
                ss[ 2] +=  2*pref *   UL[2]*UL[2];
                ss[ 3] +=    pref *(2*UL[2]*UL[2]+xiL*Alin[3*2+2]);
                ss[ 4]  =  0;
                ss[ 5]  =  0;
                ss[ 6] +=    pref *2*shear_s;
                ss[ 7] +=    pref *1*qf[3];
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
                shear_s = 2*upsilon;
                ss[ 8] +=    pref *(trA-AGA-Wg);
                ss[ 9] += -2*pref *(2*gAU+ULg*trAL-trA10);
                ss[10] +=  2*pref * U2;
                ss[11] +=    pref *(2*U2+xiL*trAL);
                ss[12]  =  0;
                ss[13]  =  0;
                ss[14] +=    pref *2*shear_s;
                ss[15] +=    pref *3*qf[3];
            }
        }
    }
    for (int j=0; j<xi.size(); ++j) {
        xi[j] *= dx;      // Convert sum to integral.
        xi[j] *= 2*M_PI;  // The azimuthal integral.
    }
    xi[0] -= 1.0;	// Calculated 1+xi, subtract 1 for xi.

    for (int j=0; j<vv.size(); ++j) {
        vv[j] *= dx;      // Convert sum to integral.
        vv[j] *= 2*M_PI;  // The azimuthal integral.
    }

    // Now sigma_perp is related to the trace by s_perp^2=(1/2)[Tr-sig_par^2]
    for (int j=8; j<ss.size(); ++j)
        ss[j] = 0.5*(ss[j] - ss[j-8]);
    for (int j=0; j<ss.size(); ++j) {
        ss[j] *= dx;	// Convert sum to integral.
        ss[j] *= 2*M_PI;	// The azimuthal integral.
    }
    std::vector<std::vector<double>> res = {xi, vv, ss};
    return res;
}

void
LSM::init(const char fname[], const double f,
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
      auto allPairs = dvsPair(rrvec[i]);
      std::vector<double> xC = allPairs[0];
      std::vector<double> vC = allPairs[1];
      std::vector<double> sC = allPairs[2];
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

double
LSM::xiRZ(const double R, const double Z, const double s2fog) {
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

std::vector<double>
LSM::xiEll(const double ss, const double s2fog,
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

void
LSM::printzFuncs(const char fbase[]) {
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

void
LSM::printqFuncs(const char fbase[]) {
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

void
LSM::printXiStuff(const char fbase[]) {
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
      std::vector<double> xC = dvsPair(rr)[0];
      fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<rr
        <<std::scientific<<std::setw(12)<<std::setprecision(4)<<qf[3];
      for (int j=0; j<xC.size(); ++j)
        fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<xC[j];
      fs<<std::endl;
    }
    fs.close();
}

void
LSM::printVpStuff(const char fbase[]) {
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
      std::vector<double> vC = dvsPair(rr)[1];
      fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<rr
        <<std::scientific<<std::setw(12)<<std::setprecision(4)<<2*qf[2];
      for (int j=0; j<vC.size(); ++j)
        fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<vC[j];
      fs<<std::endl;
    }
    fs.close();
}

void
LSM::printS2Stuff(const char fbase[]) {
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
      std::vector<double> sC = dvsPair(rr)[2];
      fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<rr;
      for (int j=0; j<sC.size(); ++j)
        fs<<std::scientific<<std::setw(12)<<std::setprecision(4)<<sC[j];
      fs<<std::endl;
    }
    fs.close();
}

