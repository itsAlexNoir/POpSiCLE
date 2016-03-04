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

## For fftw on archer
FFT_LIB	    = fftw
FFT_PATH    = ${FFTW_DIR}

MKL_LIBS   = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp

LIBS	   = ${MKL_LIBS}
OMPLIBS	   = -openmp

POPSICLE_ROOT = /home/e179/e179/alexnoir/POpSiCLE
