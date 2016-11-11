#include <cmath>

#include "gauss_legendre.hpp"

// Set the weights and abscissae, which are public variables.
void
GaussLegendre::set(const int NN) {
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
