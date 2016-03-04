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
	clean_slatec clean_debug clean_all clean \
	lib dylib check_paths

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
FFTW_OBJ        = ${addsuffix .o,${FFTW_LIST}}
SLATEC_OBJS     = ${addsuffix .o,${SLATEC_LIST}}
SHEPPACK_OBJS   = ${addsuffix .o,${SHEPPACK_LIST}}

POPSICLE_LO     = $(POPSICLE_OBJS:.o=.lo)
FFTW_LO         = $(FFTW_OBJ:.o=.lo)

ifeq (${FFT_LIB},fftw)
POPSICLE_OBJS += ${FFTW_OBJ} 
POPSICLE_LO   += ${FFTW_LO}
endif


########################################################
########################################################
##                 SET NEEDED PATHS                   ##
########################################################
########################################################

#####################
##   Set compiler ##
#####################

ifeq (${compiler},intel)
OMPLIBS = -openmp
else ifeq (${compiler},gnu)
OMPLIBS = -fopenmp
else ifeq (${compiler},cray)
OMPLIBS = 
else
OMPLIBS = 
endif

ifndef LIBTOOL
UNAME := $(shell uname)
ifeq (${UNAME},Linux)
LIBTOOL = libtool
else ifeq (${UNAME},Darwin)
LIBTOOL = glibtool
endif
endif

#####################
##   Set FFT libs  ##
#####################

ifeq (${FFT_LIB},mkl)

ifndef FFT_PATH
${error 'You must provide a path to the FFT routines.'}
endif

MKL_DFTI_OBJS = ${FFT_PATH}/include/mkl_dfti.o
MKL_DFTI_LO = ${FFT_PATH}/include/mkl_dfti.lo
ifeq (${PLIB},serial)
DIRECTIVES =-Wp,-D_USE_MKL=1,-D_USE_FFTW=0,-D_COM_MPI=0
else ifeq (${PLIB},parallel)
DIRECTIVES =-Wp,-D_USE_MKL=1,-D_USE_FFTW=0,-D_COM_MPI=1
endif

else ifeq (${FFT_LIB},fftw)

ifndef FFT_PATH
${error 'You must provide a path to the FFT routines.'}
endif

INCL += -I${FFT_PATH}/include
LIBS += -L${FFT_PATH}/lib -lfftw3 -lm
ifeq (${PLIB},serial)
DIRECTIVES =-Wp,-D_USE_MKL=0,-D_USE_FFTW=1,-D_COM_MPI=0
else ifeq (${PLIB},parallel)
DIRECTIVES =-Wp,-D_USE_MKL=0,-D_USE_FFTW=1,-D_COM_MPI=1
endif

endif

#######################################
##      Set HDF5
##  The specific libs will depend on
##  the communication choosen
########################################

INCL    += -I${HDF5_PATH}/include
LIBS    += -L${HDF5_PATH}/lib

LIBS	+= -lsz -lz -ldl -lm ${OMPLIBS}

#####################################################
#####################################################
##             Module dependencies                 ##
#####################################################
#####################################################

tools.o: constants.o
bessel.o: constants.o
gaussleg.o: constants.o
scatt_interp.o: constants.o ${SHEPPACK_OBJS}
cubic_interp.o: constants.o
scattcoords.o: constants.o scatt_interp.o
cubcoords.o: constants.o tools.o cubic_interp.o
coords.o: scattcoords.o cubcoords.o
scattboundcyl.o: constants.o gaussleg.o scatt_interp.o
scattboundcart.o: constants.o gaussleg.o scatt_interp.o
cubboundcyl.o: constants.o tools.o gaussleg.o \
	scattboundcyl.o cubic_interp.o
cubboundcart.o: constants.o tools.o gaussleg.o \
	scattboundcyl.o scattboundcart.o \
	cubic_interp.o
boundary.o: scattboundcyl.o scattboundcart.o \
	cubboundcyl.o cubboundcart.o
ifeq (${FFT_LIB},mkl)
fourier.o: ${MKL_DFTI_OBJS} constants.o tools.o
sht.o: ${MKL_DFTI_OBJS} constants.o gaussleg.o fourier.o
else ifeq (${FFT_LIB},fftw)
fourier.o: constants.o tools.o fftw.o
sht.o: constants.o gaussleg.o fftw.o fourier.o
endif
patchwork.o: constants.o tools.o
io_pop.o: constants.o patchwork.o 
io_subset.o: constants.o tools.o patchwork.o
io_surface.o: constants.o boundary.o
surface.o: constants.o tools.o boundary.o \
	io_surface.o
flux.o: constants.o gaussleg.o bessel.o sht.o \
	io_surface.o io_pop.o
observables.o: constants.o sht.o flux.o io_pop.o
popsicle.o : constants.o fourier.o \
	scattboundcyl.o scattboundcart.o \
	cubboundcyl.o  cubboundcart.o \
	scattcoords.o cubcoords.o coords.o \
	io_pop.o io_subset.o io_surface.o \
	sht.o bessel.o gaussleg.o \
	patchwork.o surface.o observables.o


########################################################

tools.lo: constants.lo
scatt_interp.lo: constants.lo ${SHEPPACK_LO}
cubic_interp.lo: constants.lo
scattcoords.lo: constants.lo scatt_interp.lo
cubcoords.lo: constants.lo tools.o cubic_interp.lo
coords.lo: scattcoords.lo cubcoords.lo
scattboundcyl.lo: constants.lo gaussleg.lo scatt_interp.lo
scattboundcart.lo: constants.lo gaussleg.lo scatt_interp.lo
cubboundcyl.lo: constants.lo tools.lo gaussleg.lo \
	scattboundcyl.lo cubic_interp.lo
cubboundcart.lo: constants.lo tools.lo gaussleg.lo \
	scattboundcyl.lo scattboundcart.lo \
	cubic_interp.lo
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
sht.lo: constants.lo gaussleg.lo fftw.lo fourier.lo
endif
patchwork.lo: constants.lo tools.lo
surface.lo: constants.lo tools.lo boundary.lo \
	io_surface.lo
flux.lo: constants.lo gaussleg.lo bessel.lo sht.lo \
	io_surface.lo io_pop.lo
observables.lo: constants.lo sht.lo flux.lo io_pop.lo
popsicle.lo : constants.lo fourier.lo \
	scattboundcyl.lo scattboundcart.lo \
	cubboundcyl.lo  cubboundcart.lo \
	scattcoords.lo cubcoords.lo coords.lo \
	io_pop.lo io_subset.lo io_surface.lo \
	sht.lo bessel.lo gaussleg.lo \
	patchwork.lo observables.lo

#################################################################
##                                                             ##
#################################################################
##                           RULES                             ##
#################################################################

lib: check_paths ${POPSICLE_OBJS}
	@echo 'Creating static library...'
	ar ${ARFLAGS} libpopsicle_${PLIB}.a ${POPSICLE_OBJS} \
	${SHEPPACK_OBJS}; \
	ranlib libpopsicle_${PLIB}.a
	cp -rf libpopsicle_${PLIB}.a ${POPSICLE_LIB}
	cp -rf *.mod ${POPSICLE_INCL}
	rm -rf libpopsicle_${PLIB}.a
	${MAKE} clean

dylib: check_paths ${POPSICLE_LO}
	@echo 'Creating dynamic library...'
	${LIBTOOL} --tag=FC --mode=link \
	${FC} ${F90FLAGS} -o libpopsicle_${PLIB}.la \
	${POPSICLE_LO} ${SHEPPACK_LO}\
	-rpath ${POPSICLE_LIB} ${LIBS}
	cp -rf *.mod ${POPSICLE_INCL}
	${MAKE} clean

install: dylib
	${LIBTOOL} --tag=FC --mode=install \
	install -c libpopsicle_${PLIB}.la \
	${POPSICLE_LIB}

check_paths:
ifndef PLIB
	${error You must specify PLIB. \
	Do you want the lib serial or parallel? }
endif
ifndef FFT_LIB
	${error 'You must provide a FFT library \
	(Do not forget the path!!'}
endif
ifneq (${FFT_LIB},mkl)
ifneq (${FFT_LIB},fftw)
	${error You must choose between mkl or \
	fftw for a FFT library.}
endif
endif
ifndef HDF5_PATH
	${error 'You must provide a path to HDF5 libraries. \
Compilation will fail miserably.'}
endif

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
	$(FC) ${F90FLAGS} -c  $< ${OMPLIBS}

%.lo:%.f90
	${LIBTOOL} --tag=FC --mode=compile \
	${FC} ${F90FLAGS} -c $< ${OMPLIBS}

## Fortran77
%.o:%.f
	${FC} ${F77FLAGS} -c $< 

%.lo:%.f
	${LIBTOOL} --tag=F77 --mode=compile \
	${FC} ${F77FLAGS} -c $<

# Pre-Processed Fortran90
%.o:%.F90
	${FC} ${F90FLAGS} ${DIRECTIVES} ${INCL} -c $<

%.lo:%.F90
	${LIBTOOL} --tag=FC --mode=compile \
	${FC} ${F90FLAGS} ${DIRECTIVES} ${INCL} -c $<
