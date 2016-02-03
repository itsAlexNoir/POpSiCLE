##*****************************************************************************
## Content:
##      Build the library
##*****************************************************************************

##-----------------------------------------------------------------
## example of use:
##  The way to compile and run all the example is simply:
##     make (or make lib, but this rule is run by default) 
##
##  In addition, you must provide the paths to some libraries
##  needed by POpSiCLE:
##  FFT_LIB: You must choose between mkl or fftw for the FFT
##  library.
##  FFT_PATH: The path to the cited library.
##  HDF5_PATH: The path to the HDF5 library. This library must be
##  the parallel version if you are building parallel POpSiCLE,
##  otherwise you need the serial one.
##  compiler: The compiler you are using. We have used the intel,
##  gnu and cray compiler to write and test POpSiCLE.
##
##  Also, you can specify more options, like:
##  FC: The compiler command (i.e. ifort)
##  F90FLAGS and F77FLAGS: Specific flags for compilation.
##
##-------------------------------------------------------------------

.PHONY: clean lib plib dylib dyplib pes_calculator help

include src/popsicle.lst

POPSICLE_LIB    = lib/libpopsicle_serial.a
POPSICLE_PLIB   = lib/libpopsicle_parallel.a
POPSICLE_DYLIB  = lib/libpopsicle_serial.la
POPSICLE_DYPLIB = lib/libpopsicle_parallel.la
POPSICLE_MOD = ${addprefix include/,${addsuffix .mod,${POPSICLE_LIST}}}


default: lib

lib:
	cd src; ${MAKE} lib "FC=${FC}" "PLIB=serial" \
		"compiler=${compiler}" "POPSICLE_ROOT=${PWD}" \
		"FFT_PATH=${FFT_PATH}" "FFT_LIB=${FFT_LIB}" \
		"F90FLAGS=${F90FLAGS}" "F77FLAGS=${F77FLAGS}" \
		"HDF5_PATH=${HDF5_PATH}"

plib:
	cd src; ${MAKE} lib "FC=${FC}" "PLIB=parallel" \
		"compiler=${compiler}" "POPSICLE_ROOT=${PWD}" \
		"FFT_PATH=${FFT_PATH}" "FFT_LIB=${FFT_LIB}" \
		"F90FLAGS=${F90FLAGS}" "F77FLAGS=${F77FLAGS}" \
		"HDF5_PATH=${HDF5_PATH}"

dylib:
	cd src; ${MAKE} dylib "FC=${FC}" "PLIB=serial" \
		"compiler=${compiler}" "POPSICLE_ROOT=${PWD}" \
		"FFT_PATH=${FFT_PATH}" "FFT_LIB=${FFT_LIB}" \
		"F90FLAGS=${F90FLAGS}" "F77FLAGS=${F77FLAGS}" \
		"HDF5_PATH=${HDF5_PATH}"

dyplib:
	cd src; ${MAKE} dylib "FC=${FC}" "PLIB=parallel" \
		"compiler=${compiler}" "POPSICLE_ROOT=${PWD}" \
		"FFT_PATH=${FFT_PATH}" "FFT_LIB=${FFT_LIB}" \
		"F90FLAGS=${F90FLAGS}" "F77FLAGS=${F77FLAGS}" \
		"HDF5_PATH=${HDF5_PATH}"

pes_calculator:
	cd utility; ${MAKE} pes_calculator "FC=${FC}" "PLIB=serial" \
		"compiler=${compiler}" "POPSICLE_ROOT=${PWD}" \
		"FFT_PATH=${FFT_PATH}" "FFT_LIB=${FFT_LIB}" \
		"F90FLAGS=${F90FLAGS}" "F77FLAGS=${F77FLAGS}" \
		"HDF5_PATH=${HDF5_PATH}"

clean:
	cd src; ${MAKE} clean
	cd examples; ${MAKE} clean
	cd bin; rm -rf *.exec
	cd lib; rm -rf *.a *.la
	cd include; rm -rf *.mod *.o

help:
	@echo
	@echo
	@echo 'Usage: make lib / make plib / make dylib / make dyplib'
	@echo '-----------------------------------------------------------------'
	@echo ' example of use:'
	@echo '  The way to compile and run all the example is simply:'
	@echo '     make (or make lib, but this rule is run by default)' 
	@echo ''
	@echo '  In addition, you must provide the paths to some libraries'
	@echo '  needed by POpSiCLE:'
	@echo '  FFT_LIB: You must choose between mkl or fftw for the FFT'
	@echo '  library.'
	@echo '  FFT_PATH: The path to the cited library.'
	@echo '  HDF5_PATH: The path to the HDF5 library. This library must be'
	@echo '  the parallel version if you are building parallel POpSiCLE,'
	@echo '  otherwise you need the serial one.'
	@echo '  compiler: The compiler you are using. We have used the intel,'
	@echo '  gnu and cray compiler to write and test POpSiCLE.'
	@echo ''
	@echo '  Also, you can specify more options, like:'
	@echo '  FC: The compiler command (i.e. ifort)'
	@echo '  F90FLAGS and F77FLAGS: Specific flags for compilation.'
	@echo ''
	@echo '-------------------------------------------------------------------'
	@echo
	@echo
