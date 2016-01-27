# -----------------------------------------------------------------------------#
#  intel15.make
# -----------------------------------------------------------------------------#


FC  = ifort
PFC = mpifort
LIBTOOL = glibtool

# Debugging flags
#OPT         = -check all
OPT         = -g -traceback -check bounds
#OPT         = -gen-interfaces -warn interfaces

# Optimisation flags
#OPT         =	-O3

# Performance flags
#OPT			= -vec-report     # For report on vectorization
#OPT			= -guide -parallel    # For report on performance
#OPT			= -guide-vec -guide-par # equivalent to the last option

F90FLAGS    =	$(OPT)
F77FLAGS    =	$(OPT)


FFT_LIB	    = mkl
FFT_PATH    = ${MKLROOT}

MKLLIBS	   = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp

LIBS	   = ${MKLLIBS}
OMPLIBS	   = -openmp

HDF5ROOT   = /opt/hdf5/serial/current
SZIPROOT   = /opt/szip/current

HDF5LIBS   = -L${HDF5ROOT}/lib/ \
                ${HDF5ROOT}/lib/libhdf5hl_fortran.a \
                ${HDF5ROOT}/lib/libhdf5_hl.a \
                ${HDF5ROOT}/lib/libhdf5_fortran.a \
                ${HDF5ROOT}/lib/libhdf5.a
HDF5LIBS	+=   -L${SZIPROOT}/lib -lsz -lz -ldl -lm
HDF5INCLUDE	=	-I${HDF5ROOT}/include
HDF5		= $(HDF5LIBS) $(HDF5INCLUDE)

PHDF5ROOT	=	/opt/hdf5/parallel/current

PHDF5LIBS    =   -L${PHDF5ROOT}/lib/ \
                ${PHDF5ROOT}/lib/libhdf5hl_fortran.a \
                ${PHDF5ROOT}/lib/libhdf5_hl.a \
                ${PHDF5ROOT}/lib/libhdf5_fortran.a \
                ${PHDF5ROOT}/lib/libhdf5.a          
PHDF5LIBS   +=   -L${SZIPROOT}/lib -lsz -lz -ldl -lm
PHDF5INCLUDE	=	-I${PHDF5ROOT}/include
PHDF5		= $(PHDF5LIBS) $(PHDF5INCLUDE)

FFTW3LIBS	+=	-I${FFTW3ROOT}/include -lfftw3 -lm

DIRECTIVES  = -Wp,-D_COM_MPI=1,-D_USE_BLAS=1

POPSICLE_ROOT = /Users/adelacalle/Documents/Codigo/POpSiCLE
