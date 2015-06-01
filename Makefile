##*****************************************************************************
## Content:
##      Build and run examples
##*****************************************************************************

.PHONY: clean


lib:
	(cd src; ${MAKE} FC="${FC}" FFLAGS="${FFLAGS}" "ARFLAGS=$(ARFLAGS)" lib)

clean:
	for d in src examples; do (cd $$d; ${MAKE} clean;); done
	/bin/rm -f *.a
