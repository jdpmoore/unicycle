# ------------------------------------------------------------------------------
# GNU gmake file.

# Change these lines to fit your system. fortranint should be 4 or 8 bytes. This
# was established when you compiled LAPACK and BLAS.

fortranint = 8
BLASLIB = -L/sw/lib -lblas
LAPACKLIB = -L/sw/lib -llapack
FORTRANLIB = -lgfortran

CPP = /sw/bin/g++
MPICPP = /sw/bin/mpicxx
FORTRAN = /sw/bin/gfortran

# Set the optimization level. I prefer '-O3'.
#opt = 
#opt = -g
opt = -O3

# Choose serial, OpenMP-parallelized, MPI-parallized versions, or hybrid
# MPI-OpenMP. I prefer OpenMP to MPI if I'm running on one shared-memory
# machine. Compression cannot be run as hybrid; only the MVP can.
#mode = s
mode = omp
#mode = mpi
#mode = hybrid

# Probably does not need to be changed:
ifeq ($(mode),s)
	MODE_FLAGS =
	ext = s
endif
ifeq ($(mode),omp)
	MODE_FLAGS = -fopenmp -DUTIL_OMP
	ext = omp
endif
ifeq ($(mode),mpi)
	CPP = $(MPICPP)
	MODE_FLAGS = -DUTIL_MPI
	ext = mpi
endif
ifeq ($(mode),hybrid)
	CPP = $(MPICPP)
	MODE_FLAGS = -fopenmp -DUTIL_OMP -DUTIL_MPI
	ext = mpi
endif

# ------------------------------------------------------------------------------
# The rest should not have to be changed.

INCLUDE = -I .
LIBS = $(LAPACKLIB) $(BLASLIB) $(FORTRANLIB)
LIBDIRS = 
OPTFLAGS = $(opt)
CPPFLAGS = $(OPTFLAGS) $(MODE_FLAGS) -DFORTRAN_INT_$(fortranint)
FFLAGS = $(OPTFLAGS) $(MODE_FLAGS) -DFORTRAN_INT_$(fortranint)
LDFLAGS = $(MODE_FLAGS)

.SUFFIXES:
.SUFFIXES: .cpp .f90 .o

CPPSRCS = src/Hd.cpp src/Compress.cpp src/Hmat.cpp src/HmatIo.cpp \
src/KeyValueFile.cpp src/CodeAnalysis.cpp src/Mpi.cpp src/CHmat.cpp \
src/SFHmat.cpp

F90SRCS = src/ambtgf.f90

OBJECTS = $(patsubst %.cpp,%.o,$(CPPSRCS)) $(patsubst %.f90,%.o,$(F90SRCS))

TGFSRC = $(patsubst %,tgf/src/%, mod_tgf.f el3dtria.f el3dtria_omp.f el3green.f el3greenh.f elh3dtria.f \
	 elh3dtria_omp.f legewhts.f mindlinpieces.f \
	 tria3adamar.f triagauc.f triahquad.f triasymq.f trilib.f)
TGFOBJ = $(patsubst %.f,%.o,$(TGFSRC)) $(patsubst %.f90,%.o,$(TGFSRC))

%.o : %.cpp
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $< -o $@

%.o : %.f
	$(FORTRAN) $(FFLAGS) -c $< -o $@

%.o : %.f90
	$(FORTRAN) $(FFLAGS) -c $< -o $@

all: libtgf libhmmvp build mvp cmvp fmvp


# Triangle Dislocation Green's Functions
libtgf: $(TGFOBJ)
	@echo "# building libtgf"
	mkdir -p lib
	ar rucs lib/libtgf_$(mode).a $(TGFOBJ)

# Library to compress and apply an H-matrix.
libhmmvp: $(OBJECTS)
	@echo "# building libhmmvp"
	mkdir -p lib
	ar rucs lib/libhmmvp_$(mode).a $(OBJECTS)

# A driver to compress an H-matrix.
build: libhmmvp libtgf src/dc3omp.o src/ambtgf.o
	mkdir -p bin
	$(CPP) src/hmmvpbuild.cpp $(INCLUDE) $(LDFLAGS) $(LIBFLAGS) $(LIBDIRS) src/dc3omp.o lib/libhmmvp_$(mode).a lib/libtgf_$(mode).a $(LIBS) -o bin/hmmvpbuild_$(mode)

# C++ examples.
mvp:
	$(CPP) examples/mvp_$(ext).cpp $(INCLUDE) $(LDFLAGS) $(LIBFLAGS) $(LIBDIRS) lib/libhmmvp_$(mode).a $(LIBS) -o examples/mvp_$(mode)

# C example.
cmvp:
ifeq ($(ext),omp)
	$(CPP) examples/cmvp_$(ext).c $(INCLUDE) $(LDFLAGS) $(LIBFLAGS) $(LIBDIRS) lib/libhmmvp_$(mode).a $(LIBS) -o examples/cmvp_$(mode)
endif

# Fortran 90 example.
fmvp:
ifeq ($(ext),omp)
	$(FORTRAN) examples/fmvp_$(ext).f90 $(OBJECTS) $(INCLUDE) $(LDFLAGS) lib/libhmmvp_$(mode).a lib/libtgf_$(mode).a $(LIBS) -lstdc++ -o examples/fmvp_$(mode)
endif

clean:
	rm -f src/*.o lib/*.a bin/* tgf/src/*.o



