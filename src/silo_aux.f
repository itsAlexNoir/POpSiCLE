CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      SUBROUTINE create_silo_file(filename)
      
      IMPLICIT NONE
      
      INCLUDE "silo.inc"
      
      INTEGER       dbfile
      INTEGER       length, ierr
      CHARACTER*(*) filename
      CHARACTER*100 name
      
      name   = trim(filename) // '.silo'
      length = len(trim(name))
c     Create the Silo file
      ierr = dbcreate(trim(name), length, DB_CLOBBER, DB_LOCAL,
     &     "A silo file created by popsicle", 31, DB_HDF5, dbfile)
      if(dbfile.eq.-1) then
         write (*, *) 'Could not create Silo file!\n'
         stop
      endif
      
      ierr = dbclose(dbfile)
      
      end subroutine create_silo_file
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      SUBROUTINE write_silo_rectilinear_mesh2D(filename, dims,
     &     meshname, x_ax, y_ax)
      
      IMPLICIT NONE
      
      INCLUDE "silo.inc"
      
      INTEGER          dbfile
      INTEGER          ierr, err
      INTEGER          dims(2), ndims, length
      DOUBLE PRECISION x_ax(dims(1))
      DOUBLE PRECISION y_ax(dims(2))
      CHARACTER*(*)    filename, meshname
      CHARACTER*100    name
      
      name   = trim(filename) // '.silo'
      length = len(trim(name))
c     Open the Silo file
      ierr = dbopen(trim(name), length, DB_HDF5, DB_APPEND, dbfile)
      if(dbfile.eq.-1) then
         write (*, *) 'Could not open Silo file!\n'
         stop
      endif
      
      ndims = 2
      length = len(trim(meshname))      
      err = dbputqm (dbfile, trim(meshname), length, "x_ax", 4, 
     &     "y_ax", 4, "z_ax", 4, x_ax, y_ax, DB_F77NULL,
     &     dims, ndims, 
     &     DB_DOUBLE, DB_COLLINEAR, DB_F77NULL, ierr)
      
      if(err.eq.-1) then
         write (*, *) 'Error writing mesh!\n'
         stop
      endif
      
      ierr = dbclose(dbfile)
      
      end subroutine write_silo_rectilinear_mesh2D

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      SUBROUTINE write_silo_rectilinear_mesh3D(filename, dims,
     &     meshname, x_ax, y_ax, z_ax )
      
      IMPLICIT NONE
      
      INCLUDE "silo.inc"
      
      INTEGER          dbfile
      INTEGER          ierr, err
      INTEGER          dims(3), ndims, length
      DOUBLE PRECISION x_ax(dims(1))
      DOUBLE PRECISION y_ax(dims(2))
      DOUBLE PRECISION z_ax(dims(3))
      CHARACTER*(*)      filename, meshname
      CHARACTER*100    name
      
      name   = trim(filename) // '.silo'
      length = len(trim(name))
c     Open the Silo file
      ierr = dbopen(trim(name), length, DB_HDF5, DB_APPEND, dbfile)
      if(dbfile.eq.-1) then
         write (*, *) 'Could not open Silo file!\n'
         stop
      endif
      
      ndims = 3
      length = len(trim(meshname))
      err = dbputqm (dbfile, trim(meshname), length, "x_ax", 4, 
     &     "y_ax", 4, "z_ax", 4, x_ax, y_ax, z_ax, dims, ndims, 
     &     DB_DOUBLE, DB_COLLINEAR, DB_F77NULL, ierr)
      
      if(err.eq.-1) then
         write (*, *) 'Error writing mesh!\n'
         stop
      endif
      
      ierr = dbclose(dbfile)
      
      END SUBROUTINE write_silo_rectilinear_mesh3D
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      SUBROUTINE write_silo_curvilinear_mesh2D(filename, dims,
     &     meshname, x_ax, y_ax )
      
      IMPLICIT NONE
      
      INCLUDE "silo.inc"
      
      INTEGER          dbfile
      INTEGER          ierr, err
      INTEGER          dims(2), ndims, length
      DOUBLE PRECISION x_ax(dims(1), dims(2))
      DOUBLE PRECISION y_ax(dims(1), dims(2))
      CHARACTER*(*)    filename, meshname
      CHARACTER*100    name
      
      
      name   = trim(filename) // '.silo'
      length = len(trim(name))
c     Open the Silo file
      ierr = dbopen(trim(name), length, DB_HDF5, DB_APPEND, dbfile)
      if(dbfile.eq.-1) then
         write (*, *) 'Could not open Silo file!\n'
         stop
      endif
      
      ndims = 2
      length = len(trim(meshname))
      err = dbputqm (dbfile, trim(meshname), length, "x_ax", 4, 
     &     "y_ax", 4, "z_ax", 4, x_ax, y_ax, DB_F77NULL,
     &     dims, ndims, DB_DOUBLE, DB_NONCOLLINEAR, DB_F77NULL, ierr)
      
      if(err.eq.-1) then
         write (*, *) 'Error writing mesh!\n'
         stop
      endif
      
      ierr = dbclose(dbfile)
      
      END SUBROUTINE write_silo_curvilinear_mesh2D

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      SUBROUTINE write_silo_curvilinear_mesh3D(filename, dims,
     &  meshname, x_ax, y_ax, z_ax )
      
      IMPLICIT NONE
      
      INCLUDE "silo.inc"
      
      INTEGER          dbfile
      INTEGER          ierr, err
      INTEGER          dims(3), ndims, length
      DOUBLE PRECISION x_ax(dims(1), dims(2), dims(3))
      DOUBLE PRECISION y_ax(dims(1), dims(2), dims(3))
      DOUBLE PRECISION z_ax(dims(1), dims(2), dims(3))
      CHARACTER*(*)    filename, meshname
      CHARACTER*100    name
      
      
      name   = trim(filename) // '.silo'
      length = len(trim(name))
c     Open the Silo file
      ierr = dbopen(trim(name), length, DB_HDF5, DB_APPEND, dbfile)
      if(dbfile.eq.-1) then
         write (*, *) 'Could not open Silo file!\n'
         stop
      endif
      
      ndims = 3
      length = len(trim(meshname))
      err = dbputqm (dbfile, trim(meshname), length, "x_ax", 4, 
     &     "y_ax", 4, "z_ax", 4, x_ax, y_ax, z_ax, dims, ndims, 
     &     DB_DOUBLE, DB_NONCOLLINEAR, DB_F77NULL, ierr)
      
      if(err.eq.-1) then
         write (*, *) 'Error writing mesh!\n'
         stop
      endif
      
      ierr = dbclose(dbfile)
      
      END SUBROUTINE write_silo_curvilinear_mesh3D
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      SUBROUTINE write_silo_function2D(filename, dims,
     &     meshname, func, funcname)
      
      IMPLICIT NONE
      
      INCLUDE "silo.inc"
      
      INTEGER          dbfile
      INTEGER          ierr, err
      INTEGER          dims(2), ndims
      INTEGER          length, lengthfunc
      DOUBLE PRECISION func(dims(1), dims(2))
      CHARACTER*(*)    filename, meshname, funcname
      CHARACTER*100    name
      
      name   = trim(filename) // '.silo'
      length = len(trim(name))
c     Open the Silo file
      ierr = dbopen(trim(name), length, DB_HDF5, DB_APPEND, dbfile)
      if(dbfile.eq.-1) then
         write (*, *) 'Could not open Silo file!\n'
         stop
      endif
      
      ndims = 2      
!     Write data
      length = len(trim(meshname))
      lengthfunc = len(trim(funcname))
      err = dbputqv1(dbfile, trim(funcname), lengthfunc, trim(meshname),
     &     length, func, dims, ndims, DB_F77NULL, 0, DB_DOUBLE,
     &     DB_NODECENT, DB_F77NULL, ierr)
      
      if(err.eq.-1) then
         write (*, *) 'Error writing 2D function!\n'
         stop
      endif
      
      ierr = dbclose(dbfile)
      
      END SUBROUTINE  write_silo_function2D
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      SUBROUTINE write_silo_function3D(filename, dims,
     &     meshname, func, funcname)
      
      IMPLICIT NONE
      
      INCLUDE "silo.inc"
      
      INTEGER          dbfile
      INTEGER          ierr, err
      INTEGER          dims(3), ndims
      INTEGER          length, lengthfunc 
      DOUBLE PRECISION func(dims(1), dims(2), dims(3))
      CHARACTER*(*)    filename, meshname, funcname
      CHARACTER*100    name
      
      name   = trim(filename) // '.silo'
      length = len(trim(name))
c     Open the Silo file
      ierr = dbopen(trim(name), length, DB_HDF5, DB_APPEND, dbfile)
      if(dbfile.eq.-1) then
         write (*, *) 'Could not open Silo file!\n'
         stop
      endif

      ndims = 3
      
!     Write data
      length = len(trim(meshname))
      lengthfunc = len(trim(funcname))
      err = dbputqv1(dbfile, trim(funcname), lengthfunc, trim(meshname),
     &     length, func, dims, ndims, DB_F77NULL, 0, DB_DOUBLE,
     &     DB_NODECENT, DB_F77NULL, ierr)
      
      if(err.eq.-1) then
         write (*, *) 'Error writing 3D function!\n'
         stop
      endif
      
      ierr = dbclose(dbfile)
      
      END SUBROUTINE  write_silo_function3D

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

