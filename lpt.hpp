#ifndef CLEFT_LPT_HPP
#define CLEFT_LPT_HPP

#include <vector>

#include "utils.hpp"

class	LPT {
public:
    void init(std::vector<double> inkLin, std::vector<double> inpLin);

    std::vector<double> Qn(const double kk) const;
    std::vector<double> Rn(const double kk) const;

private:
    std::vector<double>	kLin,pLin;
    double		dkinv;

private:
    double linInterp(const double kk) const;
    std::vector<double> Qdefs(const double r, const double x) const;
    std::vector<double> Rdefs(const double r) const;
    std::vector<double> Qnx(const double kval, const double rr) const;
};

#endif

