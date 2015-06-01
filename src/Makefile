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
POPSICLE_OBJS = constants.o coords_transformation.o \
	fourier.o interp.o popsicle.o

MKL_DFTI = ${MKLROOT}/include/mkl_dfti.o

MKL_CDFTI = ${MKLROOT}/include/mkl_cdft.o
MKL_OBJS = mkl_dfti.o mkl_cdft.o

AUX_OBJS = ${MKL_OBJS}

##-----------------------------------------------
## Module dependencies

coords_transformation.o: constants.o interp.o
fourier.o: ${MKL_DFTI} ${MKL_CDFTI} constants.o
interp.o: constants.o ${SHEPPACK_OBJS}
popsicle.o : constants.o fourier.o \
	coords_transformation.o

##-----------------------------------------------

lib: ${POPSICLE_OBJS} ${INTERP_OBJS} ${SHEPPACK_OBJS} \
	${MKL_DFTI} ${MKL_CDFTI}
	@echo 'Creating library...'
	ar ${ARFLAGS} libpopsicle.a ${POPSICLE_OBJS} \
	${INTERP_OBJS} ${SHEPPACK_OBJS}; \
	ranlib libpopsicle.a
	cp libpopsicle.a ${POPSICLE_ROOT}/lib
	cp *.mod ${SHEPPACK_DIR}/*.mod ${POPSICLE_ROOT}/lib

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
	rm -rf *.o interp/*.o ${SHEPPACK_DIR}/*.o
	rm -rf *.mod interp/*.mod ${SHEPPACK_DIR}/*.mod
	@cd ${POPSICLE_ROOT}/lib; rm -rf *.a *.mod

###############################
## Rules for makefile
###############################

## Fortran90
%.o:%.f90
	$(FC) $(F90FLAGS) -c  $<

## Fortran77
%.o:%.f
	$(FC) $(F77FLAGS) -c $<

# Pre-Processed Fortran90

%.o:%.F90
	$(FC) $(F90FLAGS) $(DIRECTIVES) -c $<
