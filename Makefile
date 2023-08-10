SDIR = source
ODIR = object
HDIR = headers
VPATH= $(SDIR)

DEPS = $(HDIR)/main.hpp $(HDIR)/ints.hpp $(HDIR)/scf.hpp $(HDIR)/powell.hpp $(HDIR)/cphf.hpp $(HDIR)/bethe.hpp

CXX = g++
CXXFLAGS = -Wall -Wcomments -O3 --std=gnu++17 -I/usr/include/eigen3 -fopenmp -lquadmath
LDFLAGS  = -fopenmp -lquadmath

_OBJ = main.o scf.o powell.o cphf.o bethe.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

xbethe: $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(ODIR)/%.o : %.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(ODIR)/*.o
	rm -f xbethe
