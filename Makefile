CXX = g++
CXXFLAGS = -Wall
CXXFLAGS += -std=c++0x
CXXFLAGS += -O2

INCLUDES = -I /usr/include -I /home/andy/armadillo/include
CXXFLAGS += ${INCLUDES}
LIBS = -L /usr/lib
LIBFLAGS = -lopenblas -llapack -lboost_system -lboost_iostreams -lboost_filesystem
DEPENDENCIES = Compound.h PotWell.h HamiltonianMaker.h MatrixSolver.h
OBJECTS = main.o Compound.o PotWell.o HamiltonianMaker.o MatrixSolver.o

all: QDots

QDots: $(OBJECTS) 
	$(CXX) $(LIBS)  $^ $(LIBFLAGS) -o  $@  

%.o: %.cxx $(DEPENDENCIES)
	$(CXX) $(CXXFLAGS) -c $^ -o $@   

clean:
	rm -f *.o

