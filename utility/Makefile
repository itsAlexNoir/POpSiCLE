##***********************************************
##
## POpSiCLE (PhOtoelectron SpeCtrum library for 
##       laser-matter intEractions)
##
##************************************************
## Content:
##    Utility makefile
##************************************************

default: pes_calculator

.PHONY: clean clean_src clean_include clean_lib \
	clean_app clean_all


ifndef POPSICLE_ROOT
POPSICLE_ROOT=..
endif
POPSICLE_SRC=${POPSICLE_ROOT}/src
POPSICLE_LIB=${POPSICLE_ROOT}/lib
POPSICLE_INCL=${POPSICLE_ROOT}/include

ifeq (${compiler},intel)
OMPLIBS = -openmp
else ifeq (${compiler},gnu)
OMPLIBS = -fopenmp
endif

HDF5_INCLUDE = ${HDF5_PATH}/include

HDF5_LIB = ${HDF5_PATH}/lib/ \
		${HDF5_PATH}/lib/libhdf5hl_fortran.a \
                ${HDF5_PATH}/lib/libhdf5_hl.a \
                ${HDF5_PATH}/lib/libhdf5_fortran.a \
                ${HDF5_PATH}/lib/libhdf5.a

ifdef SZIP_PATH
HDF5_LIB += -L${SZIP_PATH}/lib -lsz
endif


HDF5_LIB += -lz -ldl -lm

ifdef (FFT_LIB,fftw)
FFT_LIBS = -L${FFT_PATH}/lib -lfftw3 -lm
FFT_INCLUDE = -I${FFT_PATH}/include
endif

POPSICLE_LIB_SER = ${POPSICLE_LIB}/libpopsicle_serial.a
POPSICLE_LIB_PAR = ${POPSICLE_LIB}/libpopsicle_parallel.a

pes_calculator: ${POPSICLE_LIB_SER}
	${FC} ${F90FLAGS} pes_calculator.f90 -o ../bin/pes_calculator.exec \
	-I${POPSICLE_INCL} -I${HDF5_INCLUDE} ${FFT_INCLUDE} \
	${POPSICLE_LIB}/libpopsicle_serial.a -L${HDF5_LIB} \
	${FFT_LIBS} ${OMPLIBS}

${POPSICLE_LIB_SER}:
	cd ${POPSICLE_ROOT}; ${MAKE} lib "FC=${FC}" "compiler=${compiler}" \
	"FFT_LIB=${FFT_LIB}" "FFT_PATH=${FFT_PATH}" \
	"F90FLAGS=${F90FLAGS}" "F77FLAGS=${F77FLAGS}" \
	"HDF5_PATH=${HDF5ROOT}"

${POPSICLE_LIB_PAR}:
	cd ${POPSICLE_ROOT}; ${MAKE} plib "FC=${PFC}" "compiler=${compiler}" \
	"FFT_LIB=${FFT_LIB}" "FFT_PATH=${FFT_PATH}" \
	"F90FLAGS=${F90FLAGS}" "F77FLAGS=${F77FLAGS}" \
	"HDF5_PATH=${PHDF5ROOT}"

clean: clean_app

clean_all: clean_src clean_include clean_lib clean_app

clean_src:
	@echo 'Deleting source...'
	cd ${POPSICLE_SRC}; ${MAKE} clean \
	"POPSICLE_ROOT=${POPSICLE_ROOT}"

clean_lib:
	@echo 'Deleting library...'
	cd ${POPSICLE_LIB}; \
	rm -rf *.a *.la ./libs

clean_include:
	@echo 'Deleting include...'
	cd ${POPSICLE_INCL}; \
	rm -rf *.mod

clean_app:
	@echo 'Deleting app...' 
	rm -rf *.exec

###############################
## Rules for makefile
###############################

## Fortran90
%.o:%.f90
	$(FC) $(F90FLAGS) -c $<

## Fortran77
%.o:%.f
	$(FC) $(F77FLAGS) -c $<