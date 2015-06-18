##***********************************************
##
## POpSiCLE (PhOtoelectron SpeCtrum library for 
##       laser-matter intEractions)
##
##************************************************
## Content:
##    POpSiCLE examples
##************************************************

default: example

.PHONY: clean clean_src clena_lib clean_examples help

help:
	@echo "Usage: make example [function='example source name.]"
	@echo "[FC='compiler command'] [FCFLAGS='Flags for compilation']"
	@echo "[ARFLAGS='Flags for ar tool compilation']"
	@echo "[POPSICLE_ROOT='Specify path for POpSiCLE']"
	@echo "Please see example.lst file for a complete list of"
	@echo "examples available."

##-----------------------------------------------------------------
## example of use:
##  The way to compile and run all the example is simply:
##     i.e.: make example
## 	
##  If you want to compile and run a specific example
##  use function='Example name'
##      i.e.: make example function=interp1d
##
##  If you want to custom your compilation, please use
##  the compilation flags. See that all the flags have
##  default values.
## 		i.e: make example FC=ifort F90FLAGS=-03 ARFLAGS= -qv
##
##-------------------------------------------------------------------

include example.lst

ifndef ${FC}
FC := ifort
endif

ifndef threading
threading=parallel
endif

ifndef ${ARFLAGS}
64BIT=	
64BITAR=	
ARFLAGS=$(64BITAR) qv
endif

ifndef ${POPSICLE_ROOT}
POPSICLE_ROOT=..
endif
POPSICLE_SRC=${POPSICLE_ROOT}/src
POPSICLE_LIB=${POPSICLE_ROOT}/lib

ifndef ${INSTALL_DIR}
INSTALL_DIR = .
endif

ifndef ${RES_DIR}
RES_DIR = results
endif

ifndef function
function = ${EXAMPLE_LIST}
endif

ifndef ${DIRECTIVES}
DIRECTIVES=-Wp,-D_COM_MPI=0
endif

ifndef ${FFTLIB}
	FFT_LIB = fftw
	FFT_PATH = /users/adelacalle/my_libraries/fftw/current
endif

RES = $(addsuffix .res ,$(function))

POPSICLE_SER = ${POPSICLE_LIB}/libpopsicle_serial.a
POPSICLE_PAR = ${POPSICLE_LIB}/libpopsicle_parallel.a

HDF5_PATH = /users/adelacalle/my_libraries/hdf5/serial
PHDF5_PATH = /users/adelacalle/my_libraries/hdf5/parallel
SZIP_PATH = /users/adelacalle/my_libraries/szip


example: ${POPSICLE_SER} ${POPSICLE_PAR}
	${MAKE} ${RES}	

${POPSICLE_SER}:
	cd $(POPSICLE_SRC); $(MAKE) "FC=$(FC)" "F77FLAGS=$(F77FLAGS)" \
	"F90FLAGS=$(F90FLAGS)" "ARFLAGS=$(ARFLAGS)" \
	"POPSICLE_ROOT=$(POPSICLE_ROOT)" "FFT_LIB=mkl" \
	"FFT_PATH=$(MKLROOT)" "MKLROOT=$(MKLROOT)" \
	"PLIB=serial" "HDF5_PATH=$(HDF5_PATH)" \
	"SZIP_PATH=${SZIP_PATH}" "FC=ifort" \
	"DIRECTIVES=-Wp,-D_COM_MPI=0"

$(POPSICLE_PAR):
	cd $(POPSICLE_SRC); $(MAKE) "FC=$(FC)" "F77FLAGS=$(F77FLAGS)" \
	"F90FLAGS=$(F90FLAGS)" "ARFLAGS=$(ARFLAGS)" \
	"POPSICLE_ROOT=$(POPSICLE_ROOT)" "FFT_LIB=mkl" \
	"FFT_PATH=$(MKLROOT)" "MKLROOT=$(MKLROOT)" \
	"PLIB=parallel" \
	"PHDF5_PATH=$(PHDF5_PATH)" "SZIP_PATH=${SZIP_PATH}" \
	"DIRECTIVES=$(DIRECTIVES)"

${RES}: %.res: %.f90
	echo ${POPSICLE_SER}
	mkdir -p ./${RES_DIR}
	${FC}  $< -o ${RES_DIR}/$*.out ${FLAGS} -I${POPSICLE_LIB} ${POPSICLE_SER}


clean: clean_src clean_examples

clean_src:
	cd ${POPSICLE_SRC}; ${MAKE} clean \
	"POPSICLE_ROOT=${POPSICLE_ROOT}"

clean_lib:
	cd ${POPSICLE_LIB}; rm *.a *.mod

clean_examples:
	@echo 'Deleting examples...'
	@-rm -rf *.res 
	@-rm -rf *.out

###############################
## Rules for makefile
###############################

## Fortran90
%.o:%.f90
	$(FC) $(F90FLAGS) -c $<

## Fortran77
%.o:%.f
	$(FC) $(F77FLAGS) -c $<
