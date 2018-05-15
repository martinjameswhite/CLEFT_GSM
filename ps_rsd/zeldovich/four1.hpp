#ifndef	LPT_FOUR1_H_
#define	LPT_FOUR1_H_

#include	<cmath>
#include	<vector>
#include	<cstdlib>
#include	<iostream>
#include	<iomanip>
#include	"utils.hpp"




template<int N> class FFT {
// A C++ implementation of a Cooley-Tukey power-of-two complex discrete
// Fourier transform based on a routine originally by G. Rybicki.
// This is defined as a template class so that the compiler knows the
// size in advance and can optimize accordingly -- this is assumed to
// be a power of two and small enough that a highly optimized routine is
// not worthwhile.  For this reason array loops can use int's.
private:
  void fft(std::vector<double>& data, int isign) {
    // Does both the forward and inverse DFTs, depending on isign.
    // The normalization of the inverse DFT is handled below.
    int n2=2*N;
    int j=1;
    for (int i=1; i<n2; i+=2) {
      if (j > i+1) {	// Help compiler analysis of two statements below.
        double temp;
        temp=data[j-1]; data[j-1]=data[i-1]; data[i-1]=temp;
        temp=data[j+0]; data[j+0]=data[i+0]; data[i+0]=temp;
      }
      int m=N;
      while (m >= 2 && j > m) {
        j -= m;
        m >>= 1;
      }
      j += m;
    }
    for (int mmax=2; n2>mmax; mmax*=2) {
      int istep=mmax << 1;
      double theta= isign*(2*M_PI/mmax);
      double wtemp= sin(0.5*theta);
      double wpr  = -2.0*wtemp*wtemp;
      double wpi  = sin(theta);
      double wr   = 1.0;
      double wi   = 0.0;
      for (int m=1; m<mmax; m+=2) {
        for (int i=m; i<=n2; i+=istep) {
          j=i+mmax;
          double tempr=wr*data[j-1]-wi*data[j+0];
          double tempi=wr*data[j+0]+wi*data[j-1];
          data[j-1]   = data[i-1]-tempr;
          data[j+0]   = data[i+0]-tempi;
          data[i-1]  += tempr;
          data[i+0]  += tempi;
        }
        wr=(wtemp=wr)*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
      }
    }
  }
public:
  void ffft(std::vector<double>& data) {
    // Does an in-place, forward FFT.
    // Note the data array is twice the length of the FFT.
    if (data.size()!=2*N) {
      // The complex array has length twice the FFT size.
      std::cerr<<"Size mismatch in ffft."<<std::endl;
      myexit(1);
    }
    fft(data, 1);
  }
  void ifft(std::vector<double>& data) {
    // Does an in-place, inverse FFT, with normalization.
    // Note the data array is twice the length of the FFT.
    if (data.size()!=2*N) {
      // The complex array has length twice the FFT size.
      std::cerr<<"Size mismatch in ffft."<<std::endl;
      myexit(1);
    }
    fft(data,-1);
    // Apply the normalization to the inverse DFT.
    for (int nn=0; nn<data.size(); ++nn) data[nn] /= N;
  }
};
#endif
