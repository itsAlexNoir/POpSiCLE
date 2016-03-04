# -----------------------------------------------------------------------------#
#  intel15.make
# -----------------------------------------------------------------------------#


FC  = ifort
PFC = mpifort
LIBTOOL = glibtool
compiler = intel

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

# For MKL libraries
FFT_LIB	    = mkl
FFT_PATH    = ${MKLROOT}

# For FFTW libraries
#FFT_LIB	    = fftw
#FFT_PATH    = /opt/fftw/current
#FFTW_ROOT  = /opt/fftw/current
#FFTW_INCL  = -I${FFTW_ROOT}/include
#FFTW_LIBS  = -L${FFTW_ROOT}/lib -lfftw3 -lm


MKL_LIBS   = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp

LIBS	   = ${MKL_LIBS}
OMPLIBS	   = -openmp

HDF5_ROOT   = /opt/hdf5/serial/current
SZIP_ROOT   = /opt/szip/current
ZLIB_ROOT   = /opt/zlib/current

HDF5_LIBS   = -L${HDF5_ROOT}/lib/ \
                ${HDF5_ROOT}/lib/libhdf5hl_fortran.a \
                ${HDF5_ROOT}/lib/libhdf5_hl.a \
                ${HDF5_ROOT}/lib/libhdf5_fortran.a \
                ${HDF5_ROOT}/lib/libhdf5.a
HDF5_LIBS	+=   -L${SZIP_ROOT}/lib -lsz -lz -ldl -lm
HDF5_INCL  = -I${HDF5_ROOT}/include

PHDF5_ROOT   = /opt/hdf5/parallel/current

PHDF5_LIBS    =   -L${PHDF5_ROOT}/lib/ \
                ${PHDF5_ROOT}/lib/libhdf5hl_fortran.a \
                ${PHDF5_ROOT}/lib/libhdf5_hl.a \
                ${PHDF5_ROOT}/lib/libhdf5_fortran.a \
                ${PHDF5_ROOT}/lib/libhdf5.a          
PHDF5LIBS   +=   -L${SZIP_ROOT}/lib -lsz -lz -ldl -lm
PHDF5_INCL   = -I${PHDF5_ROOT}/include

POPSICLE_ROOT = /Users/adelacalle/Documents/Codigo/POpSiCLE
