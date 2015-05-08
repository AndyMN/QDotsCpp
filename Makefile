CXX=g++
CXXFLAGS=-Wall
CXXFLAGS+=-std=c++0x
CXXFLAGS+=-O2

INCLUDES=-I/usr/include -I/usr/local/include/root
LIBS=-L/usr/lib -L/opt/intel/lib/intel64 -L/opt/intel/mkl/lib/intel64 -L/usr/local/lib/root 
LIBFLAGS=-larmadillo -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 `root-config --libs`
DEPENDENCIES=Compound.h PotWell.h PotWellSolver.h
OBJECTS=main.o Compound.o PotWell.o PotWellSolver.o

all: QDots

QDots: $(OBJECTS) 
	$(CXX) $(LIBS)  $^ $(LIBFLAGS) -o  $@  

%.o: %.cxx $(DEPENDENCIES)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $^ -o $@   

clean:
	rm -f *.o

