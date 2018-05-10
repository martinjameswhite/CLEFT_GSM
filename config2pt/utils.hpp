#ifndef CLEFT_UTILS_HPP
#define CLEFT_UTILS_HPP

#include <cstdlib>
#include <iostream>
#include <exception>
#include <stdexcept>

// These can be modified later if we convert to MPI.
inline
void	myexit(const int flag) {
  exit(flag);
}

inline
void	myexception(const std::exception& e) {
  std::cout<<"Exception: "<<e.what()<<std::endl;
  myexit(1);
}

#endif

