MODULE comm_surff
  
  USE constants_pop
#if _COM_MPI
  USE MPI
#endif
  
  IMPLICIT NONE
  
  PRIVATE
  
  !---------------------------------------------------------------------------!
  !                                                                           !
  ! Public types                                                              !
  !                                                                           !
  !---------------------------------------------------------------------------!
  
  TYPE, PUBLIC :: COMMTYPE
     
     INTEGER                       :: iprocessor
     INTEGER                       :: ipk, ipt, ipp
     INTEGER                       :: numproc1dk	  
     INTEGER                       :: numproc1dt	  
     INTEGER                       :: numproc1dp	  
     INTEGER                       :: maxproc1dk	  
     INTEGER                       :: maxproc1dt	  
     INTEGER                       :: maxproc1dp	  
     INTEGER                       :: numprocessors
     INTEGER                       :: numprocessorsglobal
     INTEGER                       :: commsimulation  
     CHARACTER(LEN = 4)            :: cprocessor
     INTEGER, ALLOCATABLE          :: iparray(:, :, :)
     INTEGER, ALLOCATABLE          :: iparrayk(:)
     INTEGER, ALLOCATABLE          :: iparrayt(:)
     INTEGER, ALLOCATABLE          :: iparrayp(:)
     
  END TYPE COMMTYPE
  
  TYPE, PUBLIC :: KGRID
     
     INTEGER                       :: maxkpts, maxthetapts, maxphipts
     INTEGER                       :: numkpts, numthetapts, numphipts
     INTEGER                       :: lmax
     REAL(dp)           	    :: emin
     REAL(dp)           	    :: emax
     REAL(dp)           	    :: emax_theoretical
     REAL(dp)           	    :: deltaenergy, deltaphi
     REAL(dp)                       :: eneoffset
     REAL(dp), ALLOCATABLE	    :: energypts(:)
     REAL(dp), ALLOCATABLE	    :: kpts(:)
     REAL(dp), ALLOCATABLE	    :: thetapts(:)
     REAL(dp), ALLOCATABLE	    :: sinthetapts(:)
     REAL(dp), ALLOCATABLE	    :: costhetapts(:)
     REAL(dp), ALLOCATABLE	    :: wtheta(:)
     REAL(dp), ALLOCATABLE	    :: phipts(:)
     REAL(dp), ALLOCATABLE	    :: sinphipts(:)
     REAL(dp), ALLOCATABLE	    :: cosphipts(:)
     REAL(dp), ALLOCATABLE	    :: wphi(:)      
     
  END TYPE KGRID
  
  TYPE, PUBLIC :: RGRID
     
     INTEGER                       :: numthetapts, numphipts
     INTEGER                       :: lmax
     INTEGER                       :: numlm
     INTEGER                       :: numlmplus
     REAL(dp)           	   :: rsurface
     REAL(dp)           	   :: deltaphi
     REAL(dp), ALLOCATABLE	   :: thetapts(:)
     REAL(dp), ALLOCATABLE	   :: sinthetapts(:)
     REAL(dp), ALLOCATABLE	   :: costhetapts(:)
     REAL(dp), ALLOCATABLE         :: wtheta(:)
     REAL(dp), ALLOCATABLE         :: phipts(:)
     REAL(dp), ALLOCATABLE         :: sinphipts(:)
     REAL(dp), ALLOCATABLE         :: cosphipts(:)
     REAL(dp), ALLOCATABLE         :: wphi(:)
     LOGICAL                       :: gauge_transform
     
  END TYPE RGRID
  
  TYPE, PUBLIC :: TGRID
     
     INTEGER                       :: numtimes
     REAL(dp), ALLOCATABLE	   :: timepts(:)
     REAL(dp)                      :: deltat
     LOGICAL                       :: truncate_time_integration
     INTEGER                       :: maxtime_integrate
     LOGICAL                       :: timeaverageforrydberg
     REAL(dp)                      :: timetoaverageover
     
  END TYPE TGRID
  
  !---------------------------------------------------------------------------!
  !                                                                           !
  ! Public variables                                                          !
  !                                                                           !
  !---------------------------------------------------------------------------!
  
  TYPE(COMMTYPE), PUBLIC :: commsurff
  TYPE(KGRID), PUBLIC    :: kmesh
  TYPE(RGRID), PUBLIC    :: rmesh
  TYPE(TGRID), PUBLIC    :: tmesh
  
  !---------------------------------------------------------------------------!
  !                                                                           !
  ! Public subroutines                                                        !
  !                                                                           !
  !---------------------------------------------------------------------------!
  
  PUBLIC             :: initialize_communications
  PUBLIC             :: output_communications
  PUBLIC             :: finalize_communications
  PUBLIC             :: check_processor_count
  PUBLIC             :: initialize_ip_addresses
  PUBLIC             :: timing
  PUBLIC             :: all_processor_barrier
  PUBLIC             :: tsurff_stop
  PUBLIC             :: error
  PUBLIC             :: stopping
  PUBLIC             :: print_debug
  
  !---------------------------------------------------------------------------!
  !                                                                           !
  ! Interfaces                                                                !
  !                                                                           !
  !---------------------------------------------------------------------------!
  
  !---------------------------------------------------------------------------!
  !                                                                           !
  ! Private variables                                                         !
  !                                                                           !
  !---------------------------------------------------------------------------!
  
CONTAINS
  
  !------------------------------------------------------------------------!
  !                                                                        !
  ! Initializes the communications.                                        !
  !                                                                        !
  !------------------------------------------------------------------------!
  
  SUBROUTINE initialize_communications( numprock, numproct, numprocp )
    
    IMPLICIT NONE
    
    !--subroutine name-------------------------------------------------------!   
    
      CHARACTER(LEN = *), PARAMETER :: myname = 'initialize_communications'
      
      !--subroutine parameters-------------------------------------------------!
      
      INTEGER, INTENT(IN)           :: numprock
      INTEGER, INTENT(IN)           :: numproct
      INTEGER, INTENT(IN)           :: numprocp
      
      !--internal variables----------------------------------------------------!   
      
#if _COM_MPI
      INTEGER                       :: ip, isim, inum, myid, mysim, ierror
#endif
      
      !------------------------------------------------------------------------!
      
#if _COM_MPI
      
      ! Start up MPI
      
      CALL MPI_init( ierror )

      ! Find out number of processors.

      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, commsurff%numprocessorsglobal,       &
                          ierror )
      
      commsurff%commsimulation = MPI_COMM_WORLD
      
      ! Find out number of the processor we are working on.

      CALL MPI_COMM_RANK( commsurff%commsimulation, commsurff%iprocessor,      &
                          ierror )
          
      commsurff%numproc1dk     = numprock
      commsurff%numproc1dt     = numproct
      commsurff%numproc1dp     = numprocp
      
      commsurff%numprocessors  = commsurff%numproc1dk *		               &
                                 commsurff%numproc1dt *		               &
				 commsurff%numproc1dp
      commsurff%maxproc1dk     = commsurff%numproc1dk - 1
      commsurff%maxproc1dt     = commsurff%numproc1dt - 1
      commsurff%maxproc1dp     = commsurff%numproc1dp - 1

#else
      
      commsurff%numproc1dk     = 1
      commsurff%numproc1dt     = 1
      commsurff%numproc1dp     = 1
      
      commsurff%numprocessors  = commsurff%numproc1dk *		               &
           commsurff%numproc1dt *		                               &
           commsurff%numproc1dp
      commsurff%maxproc1dk     = commsurff%numproc1dk - 1
      commsurff%maxproc1dt     = commsurff%numproc1dt - 1
      commsurff%maxproc1dp     = commsurff%numproc1dp - 1
      
#endif

      
      WRITE(commsurff%cprocessor, '(I4.4)') commsurff%iprocessor



      END SUBROUTINE initialize_communications
        
        
      
      
      !------------------------------------------------------------------------!
      !                                                                        !
      ! Initializes the communications.                                        !
      !                                                                        !
      !------------------------------------------------------------------------!



      SUBROUTINE output_communications( )

      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER :: myname = 'output_communications'
      
      !--subroutine parameters-------------------------------------------------!
      
      !--internal variables----------------------------------------------------!   

      
      !------------------------------------------------------------------------!
      

      IF (commsurff%IPROCESSOR.EQ.0) THEN
      
         WRITE(*, *) 
      
         WRITE(*, *)                                                           &
            '--Communications Parameters-------------------------------------'
         
	 WRITE(*, *)                                     

         WRITE(*, '(A, I9)')                                                   &
            '    Number of Processors in 1D Energy       = ',                  &
	    commsurff%numproc1dk
         WRITE(*, '(A, I9)')                                                   &
            '    Number of Processors in 1D Theta        = ',                  &
	    commsurff%numproc1dt
         WRITE(*, '(A, I9)')                                                   &
            '    Number of Processors in 1D Phi          = ',                  &
	    commsurff%numproc1dp
         
         WRITE(*, *) 

	 WRITE(*, '(A, I9)')                                                   &
            '    Number of Processors                    = ',                  &
	    commsurff%numprocessors
         
         WRITE(*, *) 
      
         WRITE(*, *)                                                           &
            '----------------------------------------------------------------'

         WRITE(*, *) 
      
      ENDIF


      END SUBROUTINE output_communications
      
      
      !------------------------------------------------------------------------!
      !                                                                        !
      ! Initializes the communications.                                        !
      !                                                                        !
      !------------------------------------------------------------------------!
  
      SUBROUTINE finalize_communications( )
      
      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER :: myname = 'finalize_communications'
      
      !--subroutine parameters-------------------------------------------------!
      
      !--internal variables----------------------------------------------------!   
#if _COM_MPI
      INTEGER       :: ierror
      
      !------------------------------------------------------------------------!
      
      !     Stop MPI.      
      
      CALL MPI_finalize( ierror )
      
#else
      STOP

#endif

      END SUBROUTINE finalize_communications
      
      
      
      FUNCTION timing( )
      
      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      REAL(dp)                      :: timing
      CHARACTER(LEN = *), PARAMETER :: myname = 'timing'
      
      !--subroutine parameters-------------------------------------------------!
      
      !--internal variables----------------------------------------------------!   

      !------------------------------------------------------------------------!
      
#if _COM_MPI
      !     Stop MPI.      
      
      timing = MPI_WTIME( )
      
#endif
      
    END FUNCTION timing
    
    
    
    SUBROUTINE check_processor_count( )
      
      IMPLICIT NONE   
   
      CHARACTER(LEN = 120)  :: outputstring
      CHARACTER(LEN = 4)    :: cprocused, cproc
      

      
      IF (commsurff%numprocessorsglobal.NE.commsurff%numprocessors) THEN
         
	 
	 
	 WRITE(cprocused, '(I4.4)') commsurff%numprocessorsglobal
	 WRITE(cproc, '(I4.4)')     commsurff%numprocessors
	 
	 
	 outputstring = "Number of processors requested for job " //           &
	                "(N = " // cprocused // ") does not equal the " //     &
			"number specified in input file (N = " // cproc // ") "
	 
	 CALL tsurff_stop( outputstring )
      
      
      
      ENDIF
      
      
      
      END SUBROUTINE check_processor_count
      
      
      
      SUBROUTINE initialize_ip_addresses( )
      
      IMPLICIT NONE
      
      INTEGER       :: iprock, iproct, iprocp, ipro
      
      
      
      ALLOCATE(commsurff%iparray(0:commsurff%maxproc1dk,                       &
                                 0:commsurff%maxproc1dt,                       &
                                 0:commsurff%maxproc1dp))
      ALLOCATE(commsurff%iparrayk(0:(commsurff%numprocessors - 1)))
      ALLOCATE(commsurff%iparrayt(0:(commsurff%numprocessors - 1)))
      ALLOCATE(commsurff%iparrayp(0:(commsurff%numprocessors - 1)))

      ipro                  = -1
      commsurff%iparray  = 100.0D0
      commsurff%iparrayk = 100.0D0
      commsurff%iparrayt = 100.0D0
      commsurff%iparrayp = 100.0D0
      
      DO iprocp = 0, commsurff%maxproc1dp     
         DO iproct = 0, commsurff%maxproc1dt    
            DO iprock = 0, commsurff%maxproc1dk    
	       
	       ipro                                      = ipro + 1
	       
	       commsurff%iparray(iprock, iproct, iprocp) = ipro
	       commsurff%iparrayk(ipro)                  = iprock
	       commsurff%iparrayt(ipro)                  = iproct
	       commsurff%iparrayp(ipro)                  = iprocp
               
	       IF (commsurff%iprocessor.EQ.ipro) THEN
	       
		  commsurff%ipk = iprock
		  commsurff%ipt = iproct
		  commsurff%ipp = iprocp
	
	       ENDIF
	       
            ENDDO
         ENDDO
      ENDDO
      
      
      END SUBROUTINE initialize_ip_addresses
   
 

   !        _ _                                                   
   !   __ _| | |     _ __  _ __ ___   ___ ___  ___ ___  ___  _ __ 
   !  / _` | | |    | '_ \| '__/ _ \ / __/ _ \/ __/ __|/ _ \| '__|
   ! | (_| | | |    | |_) | | | (_) | (_|  __/\__ \__ \ (_) | |   
   !  \__,_|_|_|____| .__/|_|  \___/ \___\___||___/___/\___/|_|   
   !          |_____|_|                                           
   !  _                     _           
   ! | |__   __ _ _ __ _ __(_) ___ _ __ 
   ! | '_ \ / _` | '__| '__| |/ _ \ '__|
   ! | |_) | (_| | |  | |  | |  __/ |   
   ! |_.__/ \__,_|_|  |_|  |_|\___|_|   

      
      
      SUBROUTINE all_processor_barrier( )
      
      
      IMPLICIT NONE
   
   
   
      INTEGER       :: ierror
      
#if _COM_MPI 
      
      ! The MPI barrier routine.

      CALL MPI_barrier( commsurff%commsimulation, ierror )
     
#endif
      
      END SUBROUTINE all_processor_barrier

      
      
      !------------------------------------------------------------------------!
      !                                                                        !
      ! Error printing                                                         !
      !                                                                        !
      !------------------------------------------------------------------------!



      SUBROUTINE error( message, routine, critical )
      
      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER :: myname = 'error'
      
      !--subroutine parameters ------------------------------------------------!
      
      CHARACTER(LEN = *), INTENT(IN) :: message
      CHARACTER(LEN = *), INTENT(IN) :: routine
      LOGICAL, INTENT(IN)            :: critical
      
      !--internal variables ---------------------------------------------------!  
      
      !------------------------------------------------------------------------!



      IF (critical) THEN
	 
	 IF (commsurff%IPROCESSOR.EQ.0) THEN
	 
	    WRITE(*, *)
	    WRITE(*, *) "Critical error in subroutine: ", routine
      
         ENDIF
      
      ELSE
	 
	 IF (commsurff%IPROCESSOR.EQ.0) THEN
	 
	    WRITE(*, *)
	    WRITE(*, *) "Error message from subroutine: ", routine
       
         ENDIF
	 
      ENDIF
      
      IF (commsurff%IPROCESSOR.EQ.0) THEN
      
	 WRITE(*, *)
         WRITE(*, *) routine, ": ", message
         WRITE(*, *) 
      
      ENDIF
      
      IF (critical) THEN
	 
	 IF (commsurff%IPROCESSOR.EQ.0) THEN
	 
	    WRITE(*, *) routine, ": I will stop now."
	    WRITE(*, *)
	    WRITE(*, *)
            WRITE(*, *) 
            WRITE(*, *)                                                        &
               '----------------------------------------------------------------'
	 
	 ENDIF
	 
	 CALL finalize_communications( )
	
	 STOP
      
      ENDIF

      
      
      END SUBROUTINE error



      !------------------------------------------------------------------------!
      !                                                                        !
      ! Stop the code.                                                         !
      !                                                                        !
      !------------------------------------------------------------------------!
      
      
      
      SUBROUTINE stopping( routine, message )
      
      IMPLICIT NONE
      
      !--subroutine name-------------------------------------------------------!   
      
      CHARACTER(LEN = *), PARAMETER            :: myname = 'stopping'
      
      !--subroutine parameters ------------------------------------------------!
      
      CHARACTER(LEN = *), INTENT(IN), OPTIONAL :: message
      CHARACTER(LEN = *), INTENT(IN)           :: routine
     
      !--internal variables ---------------------------------------------------!  

      !------------------------------------------------------------------------!



      IF (commsurff%IPROCESSOR.EQ.0) THEN
	 
	 WRITE(*, *)
	 WRITE(*, *) "Stopping tsurff_calculator in subroutine: ", routine
	 
	 IF (PRESENT(message)) THEN
	 
	    WRITE(*, *)
            WRITE(*, *) routine, ": ", message
            WRITE(*, *) 
         
	 ENDIF
	 
         WRITE(*, *) 
      
         WRITE(*, *)                                                           &
            '----------------------------------------------------------------'
      ENDIF

	 
      CALL finalize_communications( )
	
      STOP

      
      
      END SUBROUTINE stopping
      
      
      
      SUBROUTINE tsurff_stop( message )
      
      IMPLICIT NONE
      
      CHARACTER(LEN = *)      :: message
      
      
      
      IF (commsurff%IPROCESSOR.EQ.0) WRITE(*, *) message
      
      CALL finalize_communications( )
      
      STOP
         
	 
  
      END SUBROUTINE tsurff_stop
      
      
      
      SUBROUTINE print_debug( message )
      
      IMPLICIT NONE
      
      CHARACTER(LEN = *)      :: message
      
      
      CALL all_processor_barrier( )
      
      IF (commsurff%iprocessor.EQ.0) WRITE(*, *) message
         
	 
  
      END SUBROUTINE print_debug
      
      
      
   END MODULE comm_surff
