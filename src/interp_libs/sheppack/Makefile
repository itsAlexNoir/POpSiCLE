#Replace with your FORTRAN 95 compiler command.
F77 = ${FC} #f95

#Place any local compiler options here.

# Use these flags if gnu compiler is used.
ifeq (${FC},'gfortran')
	FFLAGS = ${F90FLAGS} -std=f95 -pedantic -Wall -fbounds-check
else
	FFLAGS = ${F90FLAGS}
endif
#Any special compiler flags required for fixed format Fortran source.
F77FLAGS +=

VPATH = fortran95/src fortran95/drivers support
SRC = fortran95/src
DRIVERS = fortran95/drivers
SUPPORT = support

#By default, use the BLAS and/or LAPACK routines provided.
BLAS = -L. blas.o
LAPACK = -L. lapack.o
BLASDP = blas.o
LAPACKDP = lapack.o

#Uncomment the following to use your local system BLAS and/or LAPACK.
#BLAS = -lblas
#LAPACK = -llapack
#LAPACKDP = 
#BLASDP = 

#Most compilers produce lower case file names for the .mod files. 
REALMOD = real_precision.mod
LSMOD = linear_shepard_mod.mod
CSMOD = cshep2d_mod.mod
QS2MOD = qshep2d_mod.mod
QS3MOD = qshep3d_mod.mod
QSDATA = qshepmdata.mod
QSMMOD = qshepmd_mod.mod
SHMOD = sheppack.mod

#The Absoft FORTRAN compiler produces module files with capital letters
#in the file name whereas most other compilers use all lower case.
#Uncomment the following if using the Absoft compiler.
#REALMOD = REAL_PRECISION.mod
#LSMOD = LINEAR_SHEPARD_MOD.mod
#CSMOD = CSHEP2D_MOD.mod
#QS2MOD = QSHEP2D_MOD.mod
#QS3MOD = QSHEP3D_MOD.mod
#QSDATA = QSHEPMDATA.mod
#QSMMOD = QSHEPMD_MOD.mod
#SHMOD = SHEPPACK.mod


blas.o: blas.f
	$(F77) $(FFLAGS) $(F77FLAGS) -c ${SUPPORT}/blas.f

$(REALMOD): real_precision.f90
	$(F77) $(FFLAGS) -c ${SRC}/real_precision.f90

blas_shep.o: blas_shep.f90
	$(F77) $(FFLAGS) -c ${SRC}/blas_shep.f90

lapack.o: lapack.f
	$(F77) $(FFLAGS) $(F77FLAGS) -c ${SUPPORT}/lapack.f

$(LSMOD): linear_shepard.f90 $(REALMOD)
	$(F77) $(FFLAGS) -c ${SRC}/linear_shepard.f90 

$(CSMOD): cshep2d.f90 $(REALMOD)
	$(F77) $(FFLAGS) -c ${SRC}/cshep2d.f90

$(QS2MOD): qshep2d.f90 $(REALMOD)
	$(F77) $(FFLAGS) -c ${SRC}/qshep2d.f90

$(QS3MOD): qshep3d.f90 $(REALMOD)
	$(F77) $(FFLAGS) -c ${SRC}/qshep3d.f90

$(QSDATA): qshepmdata.f90 $(REALMOD)
	$(F77) $(FFLAGS) -c ${SRC}/qshepmdata.f90

$(QSMMOD): qshepmd.f90 $(REALMOD) $(QSDATA)
	$(F77) $(FFLAGS) -c ${SRC}/qshepmd.f90

$(SHMOD): sheppack.f90 $(REALMOD) $(QSMMOD) \
	$(QSDATA) $(QS3MOD) $(QS2MOD) $(CSMOD) \
	$(LSMOD)
	$(F77) $(FFLAGS) -c ${SRC}/sheppack.f90

#------------------------------------------------#

blas.lo: blas.f
	$(LIBTOOL) --tag=F77 --mode=compile \
	$(F77) $(FFLAGS) $(F77FLAGS) -c ${SUPPORT}/blas.f

$(REALMOD_DYLD): real_precision.f90
	$(LIBTOOL) --tag=FC --mode=compile \
	$(F77) $(FFLAGS) -c ${SRC}/real_precision.f90

blas_shep.lo: blas_shep.f90
	$(LIBTOOL) --tag=FC --mode=compile \
	$(F77) $(FFLAGS) -c ${SRC}/blas_shep.f90

lapack.lo: lapack.f
	$(LIBTOOL) --tag=F77 --mode=compile \
	$(F77) $(FFLAGS) $(F77FLAGS) -c ${SUPPORT}/lapack.f

$(LSMOD_DYLD): linear_shepard.f90 $(REALMOD_DYLD)
	$(LIBTOOL) --tag=FC --mode=compile \
	$(F77) $(FFLAGS) -c ${SRC}/linear_shepard.f90 

$(CSMOD_DYLD): cshep2d.f90 $(REALMOD_DYLD)
	$(LIBTOOL) --tag=FC --mode=compile \
	$(F77) $(FFLAGS) -c ${SRC}/cshep2d.f90

$(QS2MOD_DYLD): qshep2d.f90 $(REALMOD_DYLD)
	$(LIBTOOL) --tag=FC --mode=compile \
	$(F77) $(FFLAGS) -c ${SRC}/qshep2d.f90

$(QS3MOD_DYLD): qshep3d.f90 $(REALMOD_DYLD)
	$(LIBTOOL) --tag=FC --mode=compile \
	$(F77) $(FFLAGS) -c ${SRC}/qshep3d.f90

$(QSDATA_DYLD): qshepmdata.f90 $(REALMOD_DYLD)
	$(LIBTOOL) --tag=FC --mode=compile \
	$(F77) $(FFLAGS) -c ${SRC}/qshepmdata.f90

$(QSMMOD_DYLD): qshepmd.f90 $(REALMOD_DYLD) $(QSDATA_DYLD)
	$(LIBTOOL) --tag=FC --mode=compile \
	$(F77) $(FFLAGS) -c ${SRC}/qshepmd.f90

$(SHMOD_DYLD): sheppack.f90 $(REALMOD_DYLD) $(QSMMOD_DYLD) \
	$(QSDATA_DYLD) $(QS3MOD_DYLD) $(QS2MOD_DYLD) \
	$(CSMOD_DYLD) $(LSMOD_DYLD)
	$(LIBTOOL) --tag=FC --mode=compile \
	$(F77) $(FFLAGS) -c ${SRC}/sheppack.f90

#-----------------------------------------------------------#

lib: install
	ar ${ARFLAGS} libsheppack.a *.o
	ranlib libsheppack.a

dylib: install
	${LIBTOOL} --tag=FC --mode=link ${F77} \
	${FFLAGS} -o libsheppack.la *.lo

install: $(SHMOD) blas.o lapack.o blas_shep.o

install_dylib: $(SHMOD_DYLD) blas.lo lapack.lo blas_shep.lo

#---------------------------------------------------------------#

test_cshep2d: cshep2d_test.f90 $(SHMOD) blas_shep.o
	$(F77) $(FFLAGS) -o test_cshep2d ${DRIVERS}/cshep2d_test.f90 sheppack.o \
	blas_shep.o cshep2d.o

test_qshep2d: qshep2d_test.f90 $(SHMOD) blas_shep.o
	$(F77) $(FFLAGS) -o test_qshep2d ${DRIVERS}/qshep2d_test.f90 sheppack.o \
	blas_shep.o qshep2d.o

test_qshep3d: qshep3d_test.f90 $(SHMOD) blas_shep.o
	$(F77) $(FFLAGS) -o test_qshep3d ${DRIVERS}/qshep3d_test.f90 sheppack.o \
	blas_shep.o qshep3d.o

test_qshepmd: qshepmd_test.f90 $(SHMOD) blas_shep.o $(QSDATA)
	$(F77) $(FFLAGS) -o test_qshepmd ${DRIVERS}/qshepmd_test.f90 sheppack.o \
	blas_shep.o qshepmd.o ${SRC}/qshepmdata.f90

test_lshep: linear_shepard_test.f90 $(SHMOD) $(BLASDP) $(LAPACKDP)
	$(F77) $(FFLAGS) -o test_lshep ${DRIVERS}/linear_shepard_test.f90 \
	sheppack.o linear_shepard.o $(LAPACK) $(BLAS) 

test_all: test_cshep2d test_qshep2d test_qshep3d test_qshepmd test_lshep


clean: 
	rm -f test_lshep test_cshep2d test_qshep2d test_qshep3d \
	test_qshepmd *.o *.mod
	rm -rf *.lo *.la *.a ./libs
