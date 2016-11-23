#ifndef CLEFT_LSM_HPP
#define CLEFT_LSM_HPP

#include "lpt.hpp"
#include "zeldovich.hpp"
#include "spline.hpp"
#include "gauss_legendre.hpp"

// Implements the Lagrangian streaming model.
// Currently set up to take a linear power spectrum, f, F1 and F2
// at input, but could be modified to allow f, F1 and F2 to vary.
class LSM: public Zeldovich {
public:
    LSM() {}	// Do nothing.
    LSM(const char fname[], const double f,
      const double b1, const double b2, const double bs,
      const double Aeft, const double Aeft1, const double Aeft2) {
    init(fname,f,b1,b2,bs,Aeft,Aeft1,Aeft2);
    }

    void init(const char fname[], const double f,
            const double b1, const double b2, const double bs,
            const double Aeft, const double Aeft1, const double Aeft2);

    double xiRZ(const double R, const double Z, const double s2fog);

    std::vector<double> xiEll(const double ss, const double s2fog,
                            const double Apar, const double Aperp);

    void printzFuncs(const char fbase[]);
    void printqFuncs(const char fbase[]);
    void printXiStuff(const char fbase[]);
    void printVpStuff(const char fbase[]);
    void printS2Stuff(const char fbase[]);

protected:
    Spline		xispl,vvspl,stspl,spspl;
    Spline		R1spl,R2spl,Q1spl,Q2spl,Q5spl,Q8spl,Qsspl;
    std::vector<double>	X11,X22,X13,Y11,Y22,Y13,X1210,Y1210;
    std::vector<double>	V1,V3,TT;
    std::vector<double>	U1,U3,U220,U211,S2D;
    LPT			lpt;

    std::vector<double> calcEfuncs(const double q);

    void tabulateEfuncs();
    void writeSaveFile(const char fname[]);
    void readSaveFile(const char fname[]);
    std::vector<double> interpEfuncs(const double q);
    void setupQR();

    // dpair, vpair, spair.
    std::vector<std::vector<double>> dvsPair(const double rval);
};

#endif

