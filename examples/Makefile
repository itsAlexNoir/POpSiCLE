##***********************************************
##
## POpSiCLE (PhOtoelectron SpeCtrum library for 
##       laser-matter intEractions)
##
##************************************************
## Content:
##    POpSiCLE examples
##************************************************

include Sys/aigis_mkl.make

default: example_serial

.PHONY: clean clean_src clean_lib clean_examples help

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

include example_serial.lst
include example_parallel.lst

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

ifndef function_serial
function_serial = ${EXAMPLE_SER_LIST}
endif

ifndef function_parallel
function_parallel = ${EXAMPLE_PAR_LIST}
endif

FLAGS = -O3

RES_SER = $(addsuffix .res ,$(function_serial))
RES_PAR = $(addsuffix .res ,$(function_parallel))

POPSICLE_SER = ${POPSICLE_LIB}/libpopsicle_serial.a
POPSICLE_PAR = ${POPSICLE_LIB}/libpopsicle_parallel.a

example_serial: ${POPSICLE_SER}
	${MAKE} ${RES_SER}	

example_parallel: ${POPSICLE_PAR}
	${MAKE} ${RES_PAR}	

${POPSICLE_SER}:
	cd ${POPSICLE_SRC}; ${MAKE} "FC=${FC}" "F77FLAGS=${F77FLAGS}" \
	"F90FLAGS=${F90FLAGS}" "ARFLAGS=${ARFLAGS}" \
	"POPSICLE_ROOT=${POPSICLE_ROOT}" "FFT_LIB=${FFT_LIB}" \
	"FFT_PATH=${FFT_PATH}" "MKLROOT=${MKLROOT}" \
	"PLIB=serial" "HDF5_PATH=${HDF5ROOT}" \
	"SZIP_PATH=${SZIPROOT}" \
	"DIRECTIVES=-Wp,-D_COM_MPI=0"

${POPSICLE_PAR}:
	cd ${POPSICLE_SRC}; ${MAKE} "FC=${MPIFC}" "F77FLAGS=${F77FLAGS}" \
	"F90FLAGS=${F90FLAGS}" "ARFLAGS=$(ARFLAGS)" \
	"POPSICLE_ROOT=${POPSICLE_ROOT}" "FFT_LIB=${FFT_LIB}" \
	"FFT_PATH=${FFT_PATH}" "MKLROOT=${MKLROOT}" \
	"PLIB=parallel" \
	"PHDF5_PATH=${PHDF5ROOT}" "SZIP_PATH=${SZIPROOT}" \
	"DIRECTIVES=-Wp,-D_COM_MPI=1"

${RES_SER}: %.res: %.f90
	mkdir -p ./${RES_DIR}
	${FC}  $< -o ${RES_DIR}/$*.exec ${FLAGS} -I${POPSICLE_LIB} ${POPSICLE_SER} ${HDF5}

${RES_PAR}: %.res: %.f90
	mkdir -p ./${RES_DIR}
	${MPIFC}  $< -o ${RES_DIR}/$*.exec ${FLAGS} -I${POPSICLE_LIB} ${POPSICLE_PAR} ${PHDF5}


clean: clean_src clean_examples

clean_src:
	cd ${POPSICLE_SRC}; ${MAKE} clean \
	"POPSICLE_ROOT=${POPSICLE_ROOT}"

clean_lib:
	cd ${POPSICLE_LIB}; rm *.a *.mod

clean_examples:
	@echo 'Deleting examples...'
	@-rm -rf ./${RES_DIR}/*.dat 
	@-rm -rf ./${RES_DIR}/*.exec
	@-rm -rf ./${RES_DIR}/*.h5

###############################
## Rules for makefile
###############################

## Fortran90
%.o:%.f90
	$(FC) $(F90FLAGS) -c $<

## Fortran77
%.o:%.f
	$(FC) $(F77FLAGS) -c $<
