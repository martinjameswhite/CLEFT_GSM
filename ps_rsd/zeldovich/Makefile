# CXX = g++
CXX = clang++

DEPS = utils.hpp gauss_legendre.hpp spline.hpp lpt.hpp zeldovich.hpp four1.hpp

#CFLAGS += -fPIC -std=c++11 -fopenmp -g -ggdb -O3
CFLAGS += -fPIC -O3

libzeld: gauss_legendre.o spline.o lpt.o zeldovich.o wrapper.o
	${CXX} -shared -o libzeld.so gauss_legendre.o spline.o zeldovich.o wrapper.o

gauss_legendre.o: gauss_legendre.cpp ${DEPS}
	${CXX} ${CFLAGS} -c gauss_legendre.cpp

spline.o: spline.cpp ${DEPS}
	${CXX} ${CFLAGS} -c spline.cpp

zeldovich.o: zeldovich.cpp ${DEPS}
	${CXX} ${CFLAGS} -c zeldovich.cpp

wrapper.o: wrapper.cpp ${DEPS}
	${CXX} ${CFLAGS} -c wrapper.cpp

.PHONY: clean
clean:
	rm -f *.o *.so *.pyc
