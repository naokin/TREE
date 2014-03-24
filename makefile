CXX=/homec/naokin/gnu/gcc/4.8.2/bin/g++
CXXFLAGS=-g -std=c++11 -O3 -m64 -fopenmp -D_HAS_CBLAS -D_HAS_INTEL_MKL -D_DEFAULT_QUANTUM -D_ENABLE_DEFAULT_QUANTUM
#CXXFLAGS=-g -std=c++11 -O3 -m64 -fopenmp -D_HAS_CBLAS -D_HAS_INTEL_MKL -D_DEFAULT_QUANTUM -D_ENABLE_DEFAULT_QUANTUM -D_ENABLE_TTNS_DEBUG
#CXXFLAGS=-g -std=c++11 -O3 -m64 -fopenmp -D_HAS_CBLAS -D_HAS_INTEL_MKL -D_DEFAULT_QUANTUM -D_ENABLE_DEFAULT_QUANTUM -D_ENABLE_TTNS_DEBUG -D_ENABLE_BTAS_DEBUG

BLASDIR=/home100/opt/intel/mkl
BLASINCLUDE=-I$(BLASDIR)/include
BLASLIB=-L$(BLASDIR)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
#BLASLIB=-L$(BLASDIR)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lguide -lpthread

BOOSTDIR=/homec/naokin/boost/1.54.0
BOOSTINCLUDE=-I$(BOOSTDIR)/include
BOOSTLIB=-L$(BOOSTDIR)/lib -lboost_serialization

BTASDIR=/homec/naokin/btas
BTASINCLUDE=-I$(BTASDIR)/include
BTASLIB=$(BTASDIR)/lib/libbtas.a

INCLUDES=-I. $(BOOSTINCLUDE) $(BTASINCLUDE) $(BLASINCLUDE)
#LIBS=$(BLASLIB) -lgomp
LIBS=$(BOOSTLIB) $(BLASLIB)

ALLFLAGS=$(CXXFLAGS)

.C.o:
	$(CXX) $(ALLFLAGS) $(INCLUDES) -c $*.C
.f.o:
	$(F77) $(ALLFLAGS) $(INCLUDES) -c $*.f

all	: treedmrg.x



treedmrg.x	: treedmrg.o Block.o
	$(CXX) $(ALLFLAGS) $(LIBS) -o treedmrg.x treedmrg.o Block.o $(BTASLIB)

sortint.x	: sortint.o
	$(CXX) $(ALLFLAGS) $(LIBS) -o sortint.x sortint.o $(BTASLIB)

clean	:
	rm *.o; rm *.x; rm *.a

