CXX = g++
CXXFLAGS = -Wall
CXXFLAGS += -std=c++0x
CXXFLAGS += -O2

INCLUDES = -I /usr/include -I /home/andy/armadillo/include -I /home/andy/armadillo/include/armadillo_bits
CXXFLAGS += ${INCLUDES}
LIBS = -L /usr/lib -L /usr/local/lib/root 
LIBFLAGS = -lopenblas -llapack `root-config --libs`
DEPENDENCIES = Compound.h PotWell.h HamiltonianMaker.h
OBJECTS = main.o Compound.o PotWell.o HamiltonianMaker.o

all: QDots

QDots: $(OBJECTS) 
	$(CXX) $(LIBS)  $^ $(LIBFLAGS) -o  $@  

%.o: %.cxx $(DEPENDENCIES)
	$(CXX) $(CXXFLAGS) -c $^ -o $@   

clean:
	rm -f *.o

