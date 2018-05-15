#include "spline.hpp"

void
Spline::init(const std::vector<double>& x, const std::vector<double>& y) {
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

double
Spline::operator() (const double x) const {
    int	k,klo,khi;
    if (x<xa[0] || x>xa[nn-1]) {
        std::cout<<"x out of range ["<<xa[0]<<","<<xa[nn-1]<<"] in Spline." << std::endl;
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
    return(a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2[klo]+(b*b*b-b)*y2[khi])*(h*h)/6.0);
}
