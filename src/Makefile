##***********************************************
##
## POpSiCLE (PhOtoelectron SpeCtrum library for 
##       laser-matter intEractions)
##
##************************************************
## Content:
##    POpSiCLE source files
##************************************************

.PHONY: clean_src clean_lib clean_sheppack \
	clean_slatec clean_debug clean_all clean
default: lib

ifndef POPSICLE_ROOT
POPSICLE_ROOT = ..
endif

POPSICLE_LIB  = ${POPSICLE_ROOT}/lib
POPSICLE_INCL = ${POPSICLE_ROOT}/include
POPSICLE_SRC  = ${POPSICLE_ROOT}/src

SLATEC_DIR  =   interp_libs/slatec
SHEPPACK_DIR = interp_libs/sheppack
SHEPPACK_SRC = ${SHEPPACK_DIR}/fortran95/src
SHEPPACK_SPP = ${SHEPPACK_DIR}/support

VPATH = ${SLATEC_DIR} ${SHEPPACK_SRC} ${SHEPPACK_SPP}
include popsicle.lst
include ${SLATEC_DIR}/slatec.lst
include ${SHEPPACK_DIR}/sheppack.lst

POPSICLE_OBJS   = ${addsuffix .o,${POPSICLE_LIST}}
POPSICLE_LO     = ${addsuffix .lo,${POPSICLE_LIST}}
FFTW_OBJ        = ${addsuffix .o,${FFTW_LIST}}
FFTW_LO         = ${addsuffix .lo,${FFTW_LIST}}
SLATEC_OBJS     = ${addsuffix .o,${SLATEC_LIST}}
SHEPPACK_OBJS   = ${addsuffix .o,${SHEPPACK_LIST}}


ifeq (${FFT_LIB},fftw)
POPSICLE_OBJS += ${FFTW_OBJ} 
POPSICLE_LO   += ${FFTW_LO}
endif


########################################################
########################################################
##                 Set FFT libs                       ##
########################################################
########################################################

ifeq (${FFT_LIB},mkl)

ifndef FFT_PATH
${error 'You must provide a path to the FFT routines'}
endif
MKL_DFTI_OBJS = ${FFT_PATH}/include/mkl_dfti.o
MKL_DFTI_LO = ${FFT_PATH}/include/mkl_dfti.lo
ifdef DIRECTIVES
PREPROC =${DIRECTIVES},-D_USE_MKL=1,-D_USE_FFTW=0
else
DIRECTIVES=-Wp,-D_USE_MKL=1,-D_USE_FFTW=0
endif

else ifeq (${FFT_LIB},fftw)

ifndef FFT_PATH
${error You must provide a path to the FFT routines}
endif
#FFTW_OBJS = fftw.o
INCL += -I${FFT_PATH}/include
LIBS += -L${FFT_PATH}/lib -lfftw3 -lm
#F90FLAGS += -lfftw3 
ifdef DIRECTIVES
PREPROC = ${DIRECTIVES},-D_USE_MKL=0,-D_USE_FFTW=1
else
DIRECTIVES = -Wp,-D_USE_MKL=0,-D_USE_FFTW=1
endif

endif

#####################################################
##      Set libs depending on communication       ##
#####################################################

ifeq (${PLIB},parallel)
INCL += -I${PHDF5_PATH}/include
LIBS    += -L${PHDF5_PATH}/lib
else ifeq (${PLIB},serial)
ifeq (${FC},mpif90)
${error Please choose a serial compiler}
else
ifndef FC
${FC}=ifort
endif
endif	
DIRECTIVES +=,-D_COM_MPI=0
INCL = -I${HDF5_PATH}/include
LIBS    += -L${HDF5_PATH}/lib
#${error You must set PLIB variable}
endif 
#LIBS	+= -L${SZIP_PATH}/lib -lsz -lz -ldl -lm
LIBS	+= -lsz -lz -ldl -lm

#####################################################
##             Module dependencies                 ##
#####################################################

tools.o: constants.o
coords.o: constants.o scatt_interp.o
scatt_interp.o: constants.o ${SHEPPACK_OBJS}
cub_interp.o: constants.o
scattboundcyl.o: constants.o gaussleg.o scatt_interp.o
scattboundcart.o: constants.o gaussleg.o scatt_interp.o
cubboundcyl.o: constants.o tools.o gaussleg.o scattboundcyl.o \
	cubic_interp.o
cubboundcart.o: constants.o tools.o gaussleg.o scattboundcyl.o \
	scattboundcart.o cubic_interp.o
boundary.o: scattboundcyl.o scattboundcart.o\
	cubboundcyl.o cubboundcart.o
bessel.o: constants.o
gaussleg.o: constants.o
io_pop.o: constants.o patchwork.o 
io_subset.o: constants.o tools.o patchwork.o
io_surface.o: constants.o boundary.o
ifeq (${FFT_LIB},mkl)
fourier.o: ${MKL_DFTI_OBJS} constants.o tools.o
sht.o: ${MKL_DFTI_OBJS} constants.o gaussleg.o fourier.o
else ifeq (${FFT_LIB},fftw)
fourier.o: constants.o tools.o fftw.o
sht.o: constants.o gaussleg.o fftw.o
endif
patchwork.o: constants.o tools.o
surface.o: constants.o tools.o boundary.o \
	io_surface.o
flux.o: constants.o gaussleg.o bessel.o sht.o \
	io_surface.o
popsicle.o : constants.o fourier.o \
	scattboundcyl.o scattboundcart.o \
	cubboundcyl.o  cubboundcart.o \
	coords.o io_pop.o io_subset.o io_surface.o\
	sht.o bessel.o gaussleg.o patchwork.o

########################################################

tools.lo: constants.lo
coords.lo: constants.lo scatt_interp.lo
scatt_interp.lo: constants.lo ${SHEPPACK_LO}
cub_interp.lo: constants.lo
scattboundcyl.lo: constants.lo gaussleg.lo scatt_interp.lo
scattboundcart.lo: constants.lo gaussleg.lo scatt_interp.lo
cubboundcyl.lo: constants.lo tools.o gaussleg.lo \
	scattboundcyl.lo cubic_interp.lo
cubboundcart.lo: constants.lo tools.lo gaussleg.lo scattboundcyl.lo \
	scattboundcart.lo cubic_interp.lo
boundary.lo: scattboundcyl.lo scattboundcart.lo \
	cubboundcyl.lo cubboundcart.lo
bessel.lo: constants.lo
gaussleg.lo: constants.lo
io_pop.lo: constants.lo patchwork.lo 
io_subset.lo: constants.lo tools.lo patchwork.lo
io_surface.lo: constants.lo boundary.lo
ifeq (${FFT_LIB},mkl)
fourier.lo: ${MKL_DFTI_LO} constants.lo tools.lo
sht.lo: ${MKL_DFTI_LO} constants.lo gaussleg.lo fourier.lo
else ifeq (${FFT_LIB},fftw)
fourier.lo: constants.lo tools.lo fftw.lo
sht.lo: constants.lo gaussleg.lo fftw.lo
endif
patchwork.lo: constants.lo tools.lo
surface.lo: constants.lo tools.lo boundary.lo \
	io_surface.lo
flux.lo: constants.lo gaussleg.lo bessel.lo sht.lo \
	io_surface.lo
popsicle.lo : constants.lo fourier.lo \
	scattboundcyl.lo scattboundcart.lo \
	cubboundcyl.lo  cubboundcart.lo \
	coords.lo io_pop.lo io_subset.lo io_surface.lo\
	sht.lo bessel.lo gaussleg.lo patchwork.lo

#################################################################
##                                                             ##
#################################################################
##                           RULES                             ##
#################################################################

lib: ${POPSICLE_OBJS}
	@echo 'Creating static library...'
	ar ${ARFLAGS} libpopsicle_${PLIB}.a ${POPSICLE_OBJS} \
	${SHEPPACK_OBJS}; \
	ranlib libpopsicle_${PLIB}.a
	cp -rf libpopsicle_${PLIB}.a ${POPSICLE_LIB}
	cp -rf *.mod ${POPSICLE_INCL}
	rm -rf libpopsicle_${PLIB}.a

dylib: ${POPSICLE_LO}
	@echo 'Creating dynamic library...'
	${LIBTOOL} --tag=FC --mode=link \
	${FC} ${F90FLAGS} -o libpopsicle_${PLIB}.la \
	${POPSICLE_LO} ${SHEPPACK_LO}\
	-rpath ${POPSICLE_LIB} ${LIBS}
	cp -rf *.mod ${POPSICLE_INCL}

intall: dylib
	${LIBTOOL} --tag=FC --mode=install \
	install -c libpopsicle_${PLIB}.la \
	${POPSICLE_LIB}


clean: clean_debug

#clean: clean_src clean_lib
clean_all: clean_src clean_lib

clean_debug:
	rm -rf ${POPSICLE_OBJS} ${POPSICLE_LO}

clean_src:
	@echo 'Deleting source...'
	rm -rf *.a *.la *.lo *.o *.i90
	rm -rf *.mod ./libs

clean_lib:
	cd ${POPSICLE_LIB}; rm -rf *.a *.la

clean_include:
	cd ${POPSICLE_INCL}; rm -rf *.mod

clean_sheppack:
	cd ${SHEPPACK_DIR}; ${MAKE} clean

clean_slatec:
	cd ${SLATEC_DIR}; ${MAKE} clean

###############################
## Rules for makefile
###############################

## Fortran90
%.o:%.f90
	$(FC) ${F90FLAGS} -c  $<

%.lo:%.f90
	${LIBTOOL} --tag=FC --mode=compile \
	${FC} ${F90FLAGS} -c $<

## Fortran77
%.o:%.f
	${FC} ${F77FLAGS} -c $<

%.lo:%.f
	${LIBTOOL} --tag=F77 --mode=compile \
	${FC} ${F77FLAGS} -c $<

# Pre-Processed Fortran90

%.o:%.F90
	${FC} ${F90FLAGS} ${PREPROC} ${INCL} -c $<

%.lo:%.F90
	${LIBTOOL} --tag=FC --mode=compile \
	${FC} ${F90FLAGS} ${PREPROC} ${INCL} -c $<
