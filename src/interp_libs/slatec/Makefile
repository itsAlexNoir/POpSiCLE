##***********************************************
##
## SLATEC Makefile
##
##************************************************
## Content:
##    slatec source files
##************************************************

.PHONY: lib arlib clean clean_lib

include slatec.lst

SLATEC_OBJS = ${addsufix .o, ${SLATEC_LIST}}
SLATEC_LO   = ${addsuffix .lo, ${SLATEC_LIST}}

lib: ${SLATEC_OBJS}
	ar ${ARFLAGS} libslatec.a ${SLATEC_OBJS}; \
	ranlib libslatec.a

dylib: ${SLATEC_LO}
	${LIBTOOL} --tag=FC --mode=link ${FC} \
	${F90FLAGS} -o libslatec.la ${SLATEC_LO} \
	rpath ${PWD}

install: dylib
	${LIBTOOL} --tag=FC --mode=install \
	install -c libslatec.la ${PWD}

clean: clean_lib
	@echo 'Deleting .o files...'
	rm -rf *.o *.lo

clean_lib:
	@echo 'Deleting library...'
	rm -rf *.la *.a
	rm -rf ./libs

###############################
## Rules for makefile
###############################

## Fortran90
%.o:%.f90
	$(FC) ${F90FLAGS} -c  $<

%.lo:%.f90
	${LIBTOOL} --tag=FC --mode=compile ${FC} ${F90FLAGS} -c $<

## Fortran77
%.o:%.f
	${FC} ${F77FLAGS} -c $<

%.lo:%.f
	${LIBTOOL} --tag=F77 --mode=compile ${FC} ${F77FLAGS} -c $<

# Pre-Processed Fortran90

%.o:%.F90
	${FC} ${F90FLAGS} ${PREPROC} -c $< ${INCL} ${LIBS}
