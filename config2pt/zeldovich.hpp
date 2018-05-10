#ifndef CLEFT_ZELDOVICH_HPP
#define CLEFT_ZELDOVICH_HPP

#include <vector>

#include "gauss_legendre.hpp"

// Computes the correlation function within the Zel'dovich approximation.
// Based upon Carlson et al. (2013; MNRAS, 429, 1674): CLPT.
class	Zeldovich {
public:
    Zeldovich() {}
    Zeldovich(const char fname[]) { init(fname); }

    ~Zeldovich() {}

    void init(const char fname[]);
    void print_eta();
    double xiL(const double r);
    double xiZ(const double rval);
    std::vector<double> xiContributions(const double rval, const double Aeft);
    std::vector<double> v12(const double rval);
    std::vector<double> s12(const double rval);

protected:
    std::vector<double> sphBess(const double x);
    void readPowerSpectrum(const char fname[]);

    double calcSigma2();
    double calcEFTnorm();

    std::vector<double> calcQfuncs(const double q);

    double calc_nabla1(const double q);
    double calc_nabla2(const double q);

    std::vector<double> calc_Jn(const double q);

    void tabulateQfuncs();
    std::vector<double> interpQfuncs(const double q);
    std::vector<double> calcAmat(const double q[]);
    std::vector<double> calcAinv(const double q[]);

    double zeldovichIntegrand(const double r[], const double q[], const double f);

protected:
    // Table lengths for storing functions -- should be even.
    const static int	NqTable=2048,NkTable=2048;

    GaussLegendre		gl;
    std::vector<double>	kLin,pLin,etaPer,etaPar,uVal;
    std::vector<double>	xiLin,dxiLin,n2Lin,J2Lin,J3Lin,J4Lin;
    std::vector<double>	V12Lin,chiLin,zetLin;
    double		qmin,qmax,delta,dkinv,sigma2,eftNorm;
};

#endif

