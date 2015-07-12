##***********************************************
##
## POpSiCLE (PhOtoelectron SpeCtrum library for 
##       laser-matter intEractions)
##
##************************************************
## Content:
##    POpSiCLE source files
##************************************************

.PHONY: clean
default: lib

VPATH = interp interp_libs/slatec interp_libs/sheppack
SHEPPACK_DIR = interp_libs/sheppack

include interp_libs/slatec/interp_routines.lst
include interp_libs/sheppack/sheppack_routines.lst
POPSICLE_OBJS = constants.o tools.o coords_transformation.o \
	fourier.o interp.o popsicle.o io.o patchwork.o

ifeq (${FFT_LIB},fftw)
POPSICLE_OBJS += fftw.o 
endif

#################
# Set FFT libs
#################

ifeq (${FFT_LIB},mkl)
#@echo {FFT_PATH}
ifndef FFT_PATH
${error 'You must provide a path to the FFT routines'}
endif
MKL_DFTI = ${FFT_PATH}/include/mkl_dfti.o
#MKL_CDFTI = ${MKLROOT}/include/mkl_cdft.o
MKL_OBJS = mkl_dfti.o #mkl_cdft.o
ifdef DIRECTIVES
PREPROC =${DIRECTIVES},-D_USE_MKL=1,-D_USE_FFTW=0
else
DIRECTIVES=-Wp,-D_USE_MKL=1,-D_USE_FFTW=0
endif

else ifeq (${FFT_LIB},fftw)
ifndef FFT_PATH
${error You must provide a path to the FFT routines}
endif
FFTW_OBJS = fftw.o
LIBS += -L${FFT_PATH}/lib
F90FLAGS += -lfftw3 
ifdef DIRECTIVES
PREPROC = ${DIRECTIVES},-D_USE_MKL=0,-D_USE_FFTW=1
else
DIRECTIVES = -Wp,-D_USE_MKL=0,-D_USE_FFTW=1
endif

endif

########################################
# Set libs depending on communication
########################################

ifeq (${PLIB},parallel)
INCLUDS += -I${PHDF5_PATH}/include
LIBS    += -L${PHDF5_PATH}/lib #\
#	${PHDF5_PATH}/lib/libhdf5hl_fortran.a \
#	${PHDF5_PATH}/lib/libhdf5_hl.a \
#	${PHDF5_PATH}/lib/libhdf5_fortran.a \
#	${PHDF5_PATH}/lib/libhdf5.a
else ifeq (${PLIB},serial)
ifeq (${FC},mpif90)
${error Please choose a serial compiler}
else
ifndef FC
${FC}=ifort
endif
endif	
DIRECTIVES +=,-D_COM_MPI=0
INCLUDS = -I${HDF5_PATH}/include
LIBS    += -L${HDF5_PATH}/lib #\
#	${HDF5_PATH}/lib/libhdf5hl_fortran.a \
#	${HDF5_PATH}/lib/libhdf5_hl.a \
#	${HDF5_PATH}/lib/libhdf5_fortran.a \
#	${HDF5_PATH}/lib/libhdf5.a
#
#
#else
#${error You must set PLIB variable}
endif 
#LIBS	+= -L${SZIP_PATH}/lib -lsz -lz -ldl -lm
LIBS	+= -lsz -lz -ldl -lm

##-----------------------------------------------
## Module dependencies

coords_transformation.o: constants.o interp.o
interp.o: constants.o ${SHEPPACK_OBJS}
io.o: constants.o patchwork.o
popsicle.o : constants.o fourier.o \
	coords_transformation.o io.o
ifeq (${FFT_LIB},mkl)
fourier.o: ${MKL_DFTI} constants.o tools.o
else ifeq (${FFT_LIB},fftw)
fourier.o: constants.o tools.o fftw.o
endif
patchwork.o: constants.o tools.o

##-----------------------------------------------

lib: ${POPSICLE_OBJS} ${INTERP_OBJS} ${SHEPPACK_OBJS}
	@echo 'Creating library...'
	ar ${ARFLAGS} libpopsicle_${PLIB}.a ${POPSICLE_OBJS} \
	${INTERP_OBJS} ${SHEPPACK_OBJS}; \
	ranlib libpopsicle_${PLIB}.a
	cp -rf libpopsicle_${PLIB}.a ${POPSICLE_ROOT}/lib
	cp -rf *.mod ${POPSICLE_ROOT}/lib
	rm -rf libpopsicle_${PLIB}.a

${SHEPPACK_OBJS}:
	cd ${SHEPPACK_DIR}; ${MAKE} install "FC=${FC}" \
	"F90FLAGS=${F90FLAGS}";
	cp ${SHEPPACK_DIR}/*.o ${PWD}
	cp ${SHEPPACK_DIR}/*.mod ${PWD}

sheppack_tests: ${SHEPPACK_OBJS}
	cd ${SHEPPACK_DIR}; ${MAKE} test_all "FC=${FC}" \
	"F90FLAGS=${F90FLAGS}"

clean:
	@echo 'Deleting source...'
	rm -rf *.a *.o interp/*.o ${SHEPPACK_DIR}/*.o
	rm -rf *.mod interp/*.mod ${SHEPPACK_DIR}/*.mod
	@cd ${POPSICLE_ROOT}/lib; rm -rf *.a *.mod

###############################
## Rules for makefile
###############################

## Fortran90
%.o:%.f90
	$(FC) ${F90FLAGS} -c  $<

## Fortran77
%.o:%.f
	${FC} ${F77FLAGS} -c $<

# Pre-Processed Fortran90

%.o:%.F90
	${FC} ${F90FLAGS} ${PREPROC} -c $< ${INCLUDS} ${LIBS}
