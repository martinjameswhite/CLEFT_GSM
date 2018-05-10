#include	<iostream>
#include	<iomanip>

#include "utils.hpp"
#include "lsm.hpp"

int	main(int argc, char **argv)
{
  if (argc != 9) {
    std::cout<<"Usage: lesm "
             <<"<Pk-file> <ff> <b1> <b2> <bs2> <Aeft> <Aeftv> <s2FoG>"
             <<std::endl;
    myexit(1);
  }

  const double ff   = atof(argv[2]);
  const double b1   = atof(argv[3]);
  const double b2   = atof(argv[4]);
  const double bs2  = atof(argv[5]);
  const double Apar = 1;
  const double Aperp= 1;
  const double Aeft = atof(argv[6]);
  const double Aeft1= atof(argv[7]);
  const double Aeft2= 0;
  const double s2FoG= atof(argv[8]);
  try {
    LSM lsm(argv[1],ff,b1,b2,bs2,Aeft,Aeft1,Aeft2);

#ifdef	PRINTSTUFF
    lsm.printzFuncs( argv[1]);
    lsm.printqFuncs( argv[1]);
    lsm.printXiStuff(argv[1]);
    lsm.printVpStuff(argv[1]);
    lsm.printS2Stuff(argv[1]);
    return(0);
#endif

    const int Nr=61;
    const double rmin=5.0,rmax=135;
    std::vector<double> xiell;
    for (int i=0; i<Nr; ++i) {
      double rr = rmin + i*(rmax-rmin)/(Nr-1);
      xiell=lsm.xiEll(rr,s2FoG,Apar,Aperp);
      std::cout<<std::fixed<<std::setw(10)<<std::setprecision(2)<<rr
               <<std::fixed<<std::setw(15)<<std::setprecision(5)<<xiell[0]*rr*rr
               <<std::fixed<<std::setw(15)<<std::setprecision(5)<<xiell[1]*rr*rr
               <<std::endl;
    }
  } catch(std::exception& e) {myexception(e);}
  return(0);
}
