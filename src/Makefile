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

ifdef SILO_PATH
POPSICLE_OBJS += silo_aux.o
POPSICLE_LO   += silo_aux.lo
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
RANLIB	= ranlib
else ifeq (${UNAME},Darwin)
LIBTOOL = glibtool
RANLIB	= ranlib -c
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
##      Set Silo path
##  The specific libs will depend on
##  the communication choosen
########################################

ifdef SILO_PATH
DIRECTIVES +=,-D_USE_SILO=1
INCL    += -I${SILO_PATH}/include
LIBS    += -L${SILO_PATH}/lib -lsiloh5
else
DIRECTIVES +=,-D_USE_SILO=0
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
comm_surff.o: constants.o
bessel.o: constants.o
gaussleg.o: constants.o
scatt_interp.o: constants.o ${SHEPPACK_OBJS}
cubic_interp.o: constants.o
spline_interp.o: constants.o tools.o
scattcoords.o: constants.o scatt_interp.o
cubcoords.o: constants.o tools.o cubic_interp.o
coords.o: scattcoords.o cubcoords.o
waveduplicate.o: constants.o tools.o cubic_interp.o
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
else ifeq (${FFT_LIB},fftw)
fourier.o: constants.o tools.o fftw.o
endif
sharmonics.o: constants.o comm_surff.o gaussleg.o
patchwork.o: constants.o tools.o
io_pop.o: constants.o patchwork.o 
io_subset.o: constants.o tools.o patchwork.o
io_surface.o: constants.o boundary.o
surface.o: constants.o tools.o boundary.o \
	io_surface.o
corrtools.o: constants.o
corrsurface.o: constants.o tools.o corrtools.o \
	boundary.o io_surface.o io_corrsurface.o
io_corrsurface.o: constants.o boundary.o
flux.o: constants.o comm_surff.o gaussleg.o bessel.o \
	sharmonics.o tools.o io_surface.o io_pop.o
ifdef SILO_PATH
observables.o: constants.o comm_surff.o gaussleg.o \
	sharmonics.o flux.o io_pop.o silo_aux.o
else
observables.o: constants.o comm_surff.o gaussleg.o \
	sharmonics.o flux.o io_pop.o
endif
samplingpt.o: constants.o gaussleg.o fourier.o \
	io_surface.o io_pop.o
popsicle_aux.o : constants.o fourier.o comm_surff.o \
	cubic_interp.o spline_interp.o \
	scattboundcyl.o scattboundcart.o \
	cubboundcyl.o  cubboundcart.o \
	scattcoords.o cubcoords.o coords.o \
	io_pop.o io_subset.o io_surface.o \
	corrtools.o corrsurface.o io_corrsurface.o \
	sharmonics.o bessel.o gaussleg.o \
	patchwork.o surface.o samplingpt.o observables.o
popsicle.o : io_pop.o surface.o corrsurface.o \
	flux.o observables.o samplingpt.o

########################################################

tools.lo: constants.lo
comm_surff.lo: constants.lo
scatt_interp.lo: constants.lo ${SHEPPACK_LO}
cubic_interp.lo: constants.lo
spline_interp.lo: constants.lo tools.lo
scattcoords.lo: constants.lo scatt_interp.lo
cubcoords.lo: constants.lo tools.o cubic_interp.lo
coords.lo: scattcoords.lo cubcoords.lo
waveduplicate.lo: constants.lo tools.lo cubic_interp.lo
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
else ifeq (${FFT_LIB},fftw)
fourier.lo: constants.lo tools.lo fftw.lo
endif
sharmonics.lo: constants.lo comm_surff.lo gaussleg.lo
patchwork.lo: constants.lo tools.lo
surface.lo: constants.lo tools.lo boundary.lo \
	io_surface.lo
corrtools.lo: constants.lo
corrsurface.lo: constants.lo tools.lo corrtools.lo \
	boundary.lo io_surface.lo io_corrsurface.lo
io_corrsurface.lo: constants.lo boundary.lo
flux.lo: constants.lo comm_surff.lo gaussleg.lo bessel.lo \
	sharmonics.lo tools.lo io_surface.lo io_pop.lo
ifdef SILO_PATH
observables.lo: constants.lo comm_surff.lo gaussleg.lo \
	sharmonics.lo flux.lo io_pop.lo silo_aux.lo
else
observables.lo: constants.lo comm_surff.lo gaussleg.lo \
	sharmonics.lo flux.lo io_pop.lo
endif
samplingpt.lo: constants.lo gaussleg.lo fourier.lo \
	io_surface.lo io_pop.lo
popsicle_aux.lo : constants.lo fourier.lo comm_surff.lo \
	cubic_interp.lo spline_interp.lo \
	scattboundcyl.lo scattboundcart.lo \
	cubboundcyl.lo  cubboundcart.lo \
	scattcoords.lo cubcoords.lo coords.lo \
	io_pop.lo io_subset.lo io_surface.lo \
	corrtools.lo corrsurface.lo io_corrsurface.lo \
	sharmonics.lo bessel.lo gaussleg.lo \
	patchwork.lo samplingpt.lo observables.lo
popsicle.lo : io_pop.lo surface.lo corrsurface.lo \
	flux.lo observables.lo samplingpt.lo

#################################################################
##                                                             ##
#################################################################
##                           RULES                             ##
#################################################################

lib: check_paths ${POPSICLE_OBJS}
	@echo 'Creating static library...'
	ar ${ARFLAGS} libpopsicle_${PLIB}.a ${POPSICLE_OBJS} \
	${SHEPPACK_OBJS}; \
	${RANLIB} libpopsicle_${PLIB}.a
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

silo_aux.o: silo_aux.f
	${FC} ${F77FLAGS} -c $< -L${SILO_PATH}/lib \
	-I${SILO_PATH}/include -lsiloh5

# Pre-Processed Fortran90
%.o:%.F90
	${FC} ${F90FLAGS} ${DIRECTIVES} ${INCL} -c $<

%.lo:%.F90
	${LIBTOOL} --tag=FC --mode=compile \
	${FC} ${F90FLAGS} ${DIRECTIVES} ${INCL} -c $<
