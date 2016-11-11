#ifndef CLEFT_SPLINE_HPP
#define CLEFT_SPLINE_HPP

#include <vector>

#include "utils.hpp"

class	Spline {
public:
    Spline() {nn=0;}
    Spline(const std::vector<double>& x, const std::vector<double>& y) { init(x,y); }

    void init(const std::vector<double>& x, const std::vector<double>& y);
    double operator() (const double x) const;
    double xmin() const { return(xa[0]); }
    double xmax() const { return(xa[nn-1]); }

private:
    std::vector<double> xa,ya,y2;
    int	nn;
};

#endif

