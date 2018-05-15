#include	<fstream>
#include	<iomanip>
#include	<vector>

#include "utils.hpp"


int	hidden_main(const char pkfile[], const double ff, const char fname[]) {
  try {
    std::ofstream fs(fname,std::ofstream::trunc);
    if (!fs) {
      std::cerr<<"Unable to open "<<fname<<" for writing."<<std::endl;
      myexit(1);
    }
    fs<<"# Header of file."<<std::endl;
    fs<<std::fixed<<std::setw(10)<<std::setprecision(3)<<0
      <<std::fixed<<std::setw(15)<<std::setprecision(7)<<0
      <<std::fixed<<std::setw(15)<<std::setprecision(7)<<0
      <<std::endl;

    const int Nr=61;
    const double rmin=5.0,rmax=135;
    std::vector<double> xiell;
    for (int i=0; i<Nr; ++i) {
      double rr = rmin + i*(rmax-rmin)/(Nr-1);
      // xiell=lsm.xiEll(rr,s2FoG,Apar,Aperp);
      // fs<<std::fixed<<std::setw(10)<<std::setprecision(3)<<rr
      //  <<std::fixed<<std::setw(15)<<std::setprecision(7)<<xiell[0]*rr*rr
      //  <<std::fixed<<std::setw(15)<<std::setprecision(7)<<xiell[1]*rr*rr
      //  <<std::endl;
    }
    fs.close();
  } catch(std::exception& e) {myexception(e);}
  return(0);
}


extern "C" {

int	call_lesm(const char pkfile[], const double ff, const char fname[]) {
  int ret;
  ret = hidden_main(pkfile,ff,fname);
  return(ret);
}

}

