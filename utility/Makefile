##***********************************************
##
## POpSiCLE (PhOtoelectron SpeCtrum library for 
##       laser-matter intEractions)
##
##************************************************
## Content:
##    Utility makefile
##************************************************

default: tsurff_calculator

.PHONY: clean clean_src clean_include clean_lib \
	clean_app clean_all


ifndef POPSICLE_ROOT
POPSICLE_ROOT=..
endif
POPSICLE_SRC=${POPSICLE_ROOT}/src
POPSICLE_LIB=${POPSICLE_ROOT}/lib
POPSICLE_INCL=${POPSICLE_ROOT}/include

MKL_LIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp

ifeq (${compiler},intel)
OMPLIBS = -openmp
else ifeq (${compiler},gnu)
OMPLIBS = -fopenmp
endif

HDF5_INCLUDE = -I${HDF5_PATH}/include

HDF5_LIB = -L${HDF5_PATH}/lib/ \
		${HDF5_PATH}/lib/libhdf5hl_fortran.a \
                ${HDF5_PATH}/lib/libhdf5_hl.a \
                ${HDF5_PATH}/lib/libhdf5_fortran.a \
                ${HDF5_PATH}/lib/libhdf5.a

ifdef SZIP_PATH
HDF5_LIB += -L${SZIP_PATH}/lib -lsz
endif

ifdef ZLIB_PATH
HDF5_LIB += -L${ZLIB_PATH}/lib -lz
endif

HDF5_LIB += -lz -ldl -lm

ifeq (${FFT_LIB},fftw)
FFT_LIBS = -L${FFT_PATH}/lib -lfftw3 -lm
FFT_INCLUDE = -I${FFT_PATH}/include
endif

ifeq (${FFT_LIB},mkl)
FFT_LIBS = ${MKL_LIBS}
FFT_INCLUDE =
endif


POPSICLE_LIB_SER = ${POPSICLE_LIB}/libpopsicle_serial.a
POPSICLE_LIB_PAR = ${POPSICLE_LIB}/libpopsicle_parallel.a

sampling_calculator: ${POPSICLE_LIB_SER}	
	${FC} ${F90FLAGS} sampling_calculator.f90 -o ../bin/sampling_calculator.exec \
	-I${POPSICLE_INCL} ${HDF5_INCLUDE} ${FFT_INCLUDE} \
	${POPSICLE_LIB}/libpopsicle_serial.a ${HDF5_LIB} \
	${FFT_LIBS} ${OMPLIBS}

tsurff_calculator: ${POPSICLE_LIB_SER}
	${FC} ${F90FLAGS} tsurff_calculator.f90 -o ../bin/tsurff_calculator.exec \
	-I${POPSICLE_INCL} ${HDF5_INCLUDE} ${FFT_INCLUDE} \
	${POPSICLE_LIB}/libpopsicle_serial.a ${HDF5_LIB} \
	${FFT_LIBS} ${OMPLIBS}	

${POPSICLE_LIB_SER}:
	cd ${POPSICLE_ROOT}; ${MAKE} lib "FC=${FC}" "compiler=${compiler}" \
	"FFT_LIB=${FFT_LIB}" "FFT_PATH=${FFT_PATH}" \
	"F90FLAGS=${F90FLAGS}" "F77FLAGS=${F77FLAGS}" \
	"HDF5_PATH=${HDF5ROOT}" "SZIP_PATH=${SZIP_PATH}" \
	"ZLIB_PATH=${ZLIB_PATH}"

${POPSICLE_LIB_PAR}:
	cd ${POPSICLE_ROOT}; ${MAKE} plib "FC=${PFC}" "compiler=${compiler}" \
	"FFT_LIB=${FFT_LIB}" "FFT_PATH=${FFT_PATH}" \
	"F90FLAGS=${F90FLAGS}" "F77FLAGS=${F77FLAGS}" \
	"HDF5_PATH=${PHDF5ROOT}" "SZIP_PATH=${SZIP_PATH}" \
	"ZLIB_PATH=${ZLIB_PATH}"

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
