# CXX = g++
CXX = clang++

DEPS = utils.hpp gauss_legendre.hpp spline.hpp lpt.hpp zeldovich.hpp lsm.hpp

#CFLAGS = -DPRINTSTUFF

CFLAGS += -fPIC -std=c++11 -fopenmp -g -ggdb -O3

lesm: main.o libcleft
	${CXX} -fopenmp -o lesm main.o -L. -lcleft

libcleft: gauss_legendre.o spline.o lpt.o zeldovich.o lsm.o wrapper.o
	ar rcs libcleft.a gauss_legendre.o spline.o lpt.o zeldovich.o lsm.o
	${CXX} -fopenmp -shared -o libcleft.so gauss_legendre.o spline.o lpt.o zeldovich.o lsm.o wrapper.o

main.o: main.cpp ${DEPS}
	${CXX} ${CFLAGS} -c main.cpp

gauss_legendre.o: gauss_legendre.cpp ${DEPS}
	${CXX} ${CFLAGS} -c gauss_legendre.cpp

spline.o: spline.cpp ${DEPS}
	${CXX} ${CFLAGS} -c spline.cpp

lpt.o: lpt.cpp ${DEPS}
	${CXX} ${CFLAGS} -c lpt.cpp

zeldovich.o: zeldovich.cpp ${DEPS}
	${CXX} ${CFLAGS} -c zeldovich.cpp

lsm.o: lsm.cpp ${DEPS}
	${CXX} ${CFLAGS} -c lsm.cpp

wrapper.o: wrapper.cpp ${DEPS}
	${CXX} ${CFLAGS} -c wrapper.cpp

clean:
	rm *.o *.a *.so lesm
