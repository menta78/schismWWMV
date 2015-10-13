#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SAVE_BOUNDARY_INTERPOLATION_ARRAY
      USE DATAPOOL
      IMPLICIT NONE
      
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef GRIB_API_ECMWF
      SUBROUTINE INIT_GRIB_WAM_BOUNDARY
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER IFILE
      REAL(rkind) :: eTimeMjd
      LOGICAL STEPRANGE_IN
      type(FD_FORCING_GRID) :: TheInfo
      character(len=20) shortName
      integer GRIB_TYPE
      OPEN(WAV%FHNDL,FILE=WAV%FNAME,STATUS='OLD')
      WRITE(STAT%FHNDL,*) WAV%FHNDL, WAV%FNAME, BND%FHNDL, BND%FNAME
# ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_BOUND .or. (myrank .eq. 0)) THEN
# endif
        NUM_WAM_SPEC_FILES = 0
        DO
          READ( WAV%FHNDL, *, IOSTAT = ISTAT )
          IF ( ISTAT /= 0 ) EXIT
          NUM_WAM_SPEC_FILES = NUM_WAM_SPEC_FILES + 1
        END DO
        REWIND(WAV%FHNDL)
        WRITE(STAT%FHNDL,*) 'NUM_WAM_SPEC_FILES=', NUM_WAM_SPEC_FILES
        ALLOCATE(WAM_SPEC_FILE_NAMES_BND(NUM_WAM_SPEC_FILES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons_wam, allocate error 9')
        DO IFILE = 1, NUM_WAM_SPEC_FILES
          READ( WAV%FHNDL, *) WAM_SPEC_FILE_NAMES_BND(IFILE)
        END DO
        CLOSE (WAV%FHNDL)
        !
        ! reading the times
        !
        allocate(WAM_SPEC_ListTime(NUM_WAM_SPEC_FILES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_bdcons_wam, allocate error 9')
        DO IFILE = 1, NUM_WAM_SPEC_FILES
          STEPRANGE_IN = .TRUE.
          CALL READ_TIME_OF_GRIB_FILE(eTimeMjd, TRIM(WAM_SPEC_FILE_NAMES_BND(IFILE)), STEPRANGE_IN)
          WAM_SPEC_ListTime(IFILE) = eTimeMjd
        END DO
        !
        ! reading the grid
        !
        shortName='2dfd'
        GRIB_TYPE=1 ! 1 for ECMWF
        CALL READ_GRID_INFO_FROM_GRIB(TheInfo, TRIM(WAM_SPEC_FILE_NAMES_BND(1)), shortName, GRIB_TYPE)
        

        
# ifdef MPI_PARALL_GRID
      END IF
# endif
         


      
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_GRIB_WAM_BOUNDARY_WBAC(WBACOUT)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), INTENT(OUT)   :: WBACOUT(MSC,MDC,IWBMNP)
      
      END SUBROUTINE
#endif
