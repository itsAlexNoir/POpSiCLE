# -----------------------------------------------------------------------------#
#  archer_fftw.make
# -----------------------------------------------------------------------------#


FC  = ftn
PFC = ftn
LIBTOOL = libtool
compiler = cray

# Debugging flags
#OPT         = -check all
#OPT         = -g -traceback
#OPT         = -gen-interfaces -warn interfaces

# Optimisation flags
OPT         =	-O3

# Performance flags
#OPT			= -vec-report     # For report on vectorization
#OPT			= -guide -parallel    # For report on performance
#OPT			= -guide-vec -guide-par # equivalent to the last option

F90FLAGS    =	$(OPT)
F77FLAGS    =	$(OPT)


FFT_LIB	    = fftw
FFT_PATH    = ${FFTW_DIR}

MKLLIBS	   = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp

LIBS	   = ${MKLLIBS}
OMPLIBS	   = -openmp

SZIPROOT   = /users/adelacalle/my_libraries/szip

HDF5LIBS   = -L${HDF5ROOT}/lib/ \
                ${HDF5ROOT}/lib/libhdf5hl_fortran.a \
                ${HDF5ROOT}/lib/libhdf5_hl.a \
                ${HDF5ROOT}/lib/libhdf5_fortran.a \
                ${HDF5ROOT}/lib/libhdf5.a
HDF5LIBS	+=   -L${SZIPROOT}/lib -lsz -lz -ldl -lm
HDF5INCLUDE	=	-I${HDF5ROOT}/include
HDF5         = $(HDF5LIBS) $(HDF5INCLUDE)

PHDF5ROOT    = ${HDF5_ROOT}

PHDF5LIBS    = -L${PHDF5ROOT}/lib/ \
                ${PHDF5ROOT}/lib/libhdf5hl_fortran.a \
                ${PHDF5ROOT}/lib/libhdf5_hl.a \
                ${PHDF5ROOT}/lib/libhdf5_fortran.a \
                ${PHDF5ROOT}/lib/libhdf5.a          
PHDF5LIBS   += -L${SZIPROOT}/lib -lsz -lz -ldl -lm
PHDF5INCLUDE = -I${PHDF5ROOT}/include
PHDF5     = $(PHDF5LIBS) $(PHDF5INCLUDE)

FFTW3LIBS += -I${FFTW3ROOT}/include -lfftw3 -lm

POPSICLE_ROOT = /home/e179/e179/alexnoir/POpSiCLE
