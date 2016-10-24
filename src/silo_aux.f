CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      SUBROUTINE create_silo_file(filename)
      
      IMPLICIT NONE
      
      INCLUDE "silo.inc"

      INTEGER  dbfile
      INTEGER  length, ierr
      CHARACTER*100 filename, name

      name   = trim(filename) // '.silo'
      length = len(name)
c     Create the Silo file
      ierr = dbcreate(name, length, DB_CLOBBER, DB_LOCAL,
     &     "A silo file created by popsicle", 31, DB_HDF5, dbfile)
      if(dbfile.eq.-1) then
         write (*, *) 'Could not create Silo file!\n'
         stop
      endif
      
      ierr = dbclose(dbfile)

      end subroutine create_silo_file
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      SUBROUTINE write_silo_rectilinear_mesh2D(filename, dims,
     &     x_ax, y_ax)
      
      IMPLICIT NONE
      
      INCLUDE "silo.inc"

      INTEGER          dbfile
      INTEGER          ierr, err
      INTEGER          dims(:), ndims, length
      DOUBLE PRECISION x_ax(:)
      DOUBLE PRECISION y_ax(:)
      CHARACTER*100 filename, name
      
      
      name   = trim(filename) // '.silo'
      length = len(name)
c     Open the Silo file
      ierr = dbopen(name, length, DB_HDF5, DB_APPEND, dbfile)
      if(dbfile.eq.-1) then
         write (*, *) 'Could not open Silo file!\n'
         stop
      endif

      ndims = 2
      err = dbputqm (dbfile, "curvimesh3D", 11, "x_ax", 4, 
     &     "y_ax", 4, x_ax, y_ax, dims, ndims, 
     &     DB_FLOAT, DB_COLLINEAR, DB_F77NULL, ierr)

      if(err.eq.-1) then
         write (*, *) 'Error writing mesh!\n'
         stop
      endif
            
      ierr = dbclose(dbfile)
      
      end subroutine write_silo_rectilinear_mesh2D

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      SUBROUTINE write_silo_rectilinear_mesh3D(filename, dims,
     &     x_ax, y_ax, z_ax )
      
      IMPLICIT NONE
      
      INCLUDE "silo.inc"
      
      INTEGER          dbfile
      INTEGER          ierr, err
      INTEGER          dims(:), ndims, length
      DOUBLE PRECISION x_ax(:)
      DOUBLE PRECISION y_ax(:)
      DOUBLE PRECISION z_ax(:)
      CHARACTER*100 filename, name
      
       
      name   = trim(filename) // '.silo'
      length = len(name)
c     Open the Silo file
      ierr = dbopen(name, length, DB_HDF5, DB_APPEND, dbfile)
      if(dbfile.eq.-1) then
         write (*, *) 'Could not open Silo file!\n'
         stop
      endif

      ndims = 3
      err = dbputqm (dbfile, "curvimesh3D", 11, "x_ax", 4, 
     &     "y_ax", 4, "z_ax", 4, x_ax, y_ax, z_ax, dims, ndims, 
     &     DB_FLOAT, DB_COLLINEAR, DB_F77NULL, ierr)
      
      if(err.eq.-1) then
         write (*, *) 'Error writing mesh!\n'
         stop
      endif
      
      ierr = dbclose(dbfile)

      END SUBROUTINE write_silo_rectilinear_mesh3D
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      SUBROUTINE write_silo_curvilinear_mesh2D(filename, dims,
     &     x_ax, y_ax )
      
      IMPLICIT NONE
      
      INCLUDE "silo.inc"
      
      INTEGER          dbfile
      INTEGER          ierr, err
      INTEGER          dims(:), ndims, length
      DOUBLE PRECISION x_ax(:, :)
      DOUBLE PRECISION y_ax(:, :)
      CHARACTER*100 filename, name
      
      
      name   = trim(filename) // '.silo'
      length = len(name)
c     Open the Silo file
      ierr = dbopen(name, length, DB_HDF5, DB_APPEND, dbfile)
      if(dbfile.eq.-1) then
         write (*, *) 'Could not open Silo file!\n'
         stop
      endif

      ndims = 2
      err = dbputqm (dbfile, "curvimesh2D", 11, "x_ax", 4, 
     &     "y_ax", 4, x_ax, y_ax, dims, ndims, 
     &     DB_FLOAT, DB_NONCOLLINEAR, DB_F77NULL, ierr)

      if(err.eq.-1) then
         write (*, *) 'Error writing mesh!\n'
         stop
      endif
      
      ierr = dbclose(dbfile)

      END SUBROUTINE write_silo_curvilinear_mesh2D

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      SUBROUTINE write_silo_curvilinear_mesh3D(filename, dims,
     &  x_ax, y_ax, z_ax )
      
      IMPLICIT NONE
      
      INCLUDE "silo.inc"
      
      INTEGER          dbfile
      INTEGER          ierr, err
      INTEGER          dims(:), ndims, length
      DOUBLE PRECISION x_ax(:, :, :)
      DOUBLE PRECISION y_ax(:, :, :)
      DOUBLE PRECISION z_ax(:, :, :)
      CHARACTER*100 filename, name
      
      
      name   = trim(filename) // '.silo'
      length = len(name)
c     Open the Silo file
      ierr = dbopen(name, length, DB_HDF5, DB_APPEND, dbfile)
      if(dbfile.eq.-1) then
         write (*, *) 'Could not open Silo file!\n'
         stop
      endif

      ndims = 3
      err = dbputqm (dbfile, "curvimesh3D", 11, "x_ax", 4, 
     &     "y_ax", 4, "z_ax", 4, x_ax, y_ax, z_ax, dims, ndims, 
     &     DB_FLOAT, DB_NONCOLLINEAR, DB_F77NULL, ierr)
      
      if(err.eq.-1) then
         write (*, *) 'Error writing mesh!\n'
         stop
      endif
      
      ierr = dbclose(dbfile)
      
      END SUBROUTINE write_silo_curvilinear_mesh3D
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      SUBROUTINE write_silo_function2D(filename, dims,
     &     func)
      
      IMPLICIT NONE
      
      INCLUDE "silo.inc"
      
      INTEGER          dbfile
      INTEGER          ierr, err
      INTEGER          dims(:), ndims, length
      DOUBLE PRECISION func(:, :)
      CHARACTER*100 filename, name
      
      
      name   = trim(filename) // '.silo'
      length = len(name)
c     Open the Silo file
      ierr = dbopen(name, length, DB_HDF5, DB_APPEND, dbfile)
      if(dbfile.eq.-1) then
         write (*, *) 'Could not open Silo file!\n'
         stop
      endif

      ndims = 2
      
!     Write data
      err = dbputqv1(dbfile, "function2D", 10, "mesh", 4, func, dims, 
     &     ndims, DB_F77NULL, 0, DB_FLOAT, DB_NODECENT, DB_F77NULL,
     &     ierr)
      
      if(err.eq.-1) then
         write (*, *) 'Error writing mesh!\n'
         stop
      endif
      
      ierr = dbclose(dbfile)
      
      END SUBROUTINE  write_silo_function2D
      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      SUBROUTINE write_silo_function3D(filename, dims,
     &     func)
      
      IMPLICIT NONE
      
      INCLUDE "silo.inc"
      
      INTEGER          dbfile
      INTEGER          ierr, err
      INTEGER          dims(:), ndims, length 
      DOUBLE PRECISION func(:, :, :)
      CHARACTER*100 filename, name
      
      
      name   = trim(filename) // '.silo'
      length = len(name)
c     Open the Silo file
      ierr = dbopen(name, length, DB_HDF5, DB_APPEND, dbfile)
      if(dbfile.eq.-1) then
         write (*, *) 'Could not open Silo file!\n'
         stop
      endif

      ndims = 3
      
!     Write data
      err = dbputqv1(dbfile, "function3D", 10, "mesh", 4, func, dims, 
     &     ndims, DB_F77NULL, 0, DB_FLOAT, DB_NODECENT, DB_F77NULL,
     &     ierr)
      
      if(err.eq.-1) then
         write (*, *) 'Error writing mesh!\n'
         stop
      endif
      
      ierr = dbclose(dbfile)
      
      END SUBROUTINE  write_silo_function3D

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
