# -----------------------------------------------------------------------------#
#  lobster_fftw.make
# -----------------------------------------------------------------------------#


FC  = ifort
MPIF90 = mpif90

# Debugging flags
#OPT         = -check all
#OPT         = -g -traceback
#OPT         = -gen-interfaces -warn interfaces

# Optimisation flags
OPT         =	-O3 -mkl=sequential # -heap-arrays

# Performance flags
#OPT			= -vec-report     # For report on vectorization
#OPT			= -guide -parallel    # For report on performance
#OPT			= -guide-vec -guide-par # equivalent to the last option

F90FLAGS    =	$(OPT)
F77FLAGS    =	$(OPT)


FFT_LIB	    = fftw
FFT_PATH    = /usr/local/codes/fftw/3.3.4/intel/14.0

#MKLROOT    = ${MKLROOT}

MKLLIBS	   = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp

LIBS	   = ${MKLLIBS}
OMPLIBS	   = -openmp

#HDF5ROOT   = /users/adelacalle/my_libraries/hdf5/serial
SZIPROOT   = /s/adelacalle/my_libs/szip/current

HDF5LIBS   = -L${HDF5ROOT}/lib/ \
                ${HDF5ROOT}/lib/libhdf5hl_fortran.a \
                ${HDF5ROOT}/lib/libhdf5_hl.a \
                ${HDF5ROOT}/lib/libhdf5_fortran.a \
                ${HDF5ROOT}/lib/libhdf5.a
HDF5LIBS	+=   -L${SZIPROOT}/lib -lsz -lz -ldl -lm
HDF5INCLUDE	=	-I${HDF5ROOT}/include
HDF5		= $(HDF5LIBS) $(HDF5INCLUDE)

PHDF5ROOT	= /s/adelacalle/my_libs/hdf5/current

PHDF5LIBS    =   -L${PHDF5ROOT}/lib/ \
                ${PHDF5ROOT}/lib/libhdf5hl_fortran.a \
                ${PHDF5ROOT}/lib/libhdf5_hl.a \
                ${PHDF5ROOT}/lib/libhdf5_fortran.a \
                ${PHDF5ROOT}/lib/libhdf5.a          
PHDF5LIBS   +=   -L${SZIPROOT}/lib -lsz -lz -ldl -lm
PHDF5INCLUDE	=	-I${PHDF5ROOT}/include
PHDF5		= $(PHDF5LIBS) $(PHDF5INCLUDE)

FFTW3ROOT	= /usr/local/codes/fftw/3.3.4/intel/14.0
#FFTW3LIBS	=	-L${FFTW3ROOT}/lib
FFTW3LIBS	+=	-I${FFTW3ROOT}/include -lfftw3 -lm

DIRECTIVES  = -Wp,-D_COM_MPI=1,-D_USE_BLAS=1
