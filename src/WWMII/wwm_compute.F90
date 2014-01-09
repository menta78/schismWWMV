#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_SIMPLE_EXPLICIT
        USE DATAPOOL
        IMPLICIT NONE

        INTEGER           :: IS, ID, IP
        REAL(rkind)       :: VEC2RAD, DEG

#ifdef TIMINGS
        REAL(rkind)       :: TIME1, TIME2, TIME3, TIME4, TIME5
        REAL(rkind)       :: TIME6, TIME7, TIME8, TIME9, TIME10, TIME11, TIME12, TIME13
#endif

         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'START COMPUTE COMPUTE_SIMPLE_EXPLICIT'
         FLUSH(STAT%FHNDL)

         AC1 = AC2 
         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER ENTERING COMPUTE ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN COMPUTE 1') 
         ENDIF

         IF (.NOT. LSTEA .AND. .NOT. LQSTEA) THEN
           DT4A = MAIN%DELT
           DT4S = DT4A
           DT4D = 0.5_rkind*DT4A
           DT4F = 0.5_rkind*DT4A 
         ELSE IF (LQSTEA) THEN
           DT4A = DT_ITER
           DT4S = DT4A
           DT4D = 0.5_rkind*DT4A
           DT4F = 0.5_rkind*DT4A
         END IF

#ifdef TIMINGS
         CALL MY_WTIME(TIME1)
#endif

         CALL COMPUTE_DIFFRACTION

#ifdef TIMINGS
         CALL MY_WTIME(TIME2)
#endif

         IF (FMETHOD .GT. 0) CALL COMPUTE_FREQUENCY()
         IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION()

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER DIRECTION AND FREQUENCY -1- ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN COMPUTE 2')
         ENDIF
  
#ifdef TIMINGS
         CALL MY_WTIME(TIME3)
#endif
         IF (FMETHOD .GT. 0) CALL COMPUTE_FREQUENCY()
         IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION()

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER DIRECTION AND FREQUENCY -2- ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN COMPUTE 3') 
         ENDIF

#ifdef TIMINGS
         CALL MY_WTIME(TIME4)
#endif
         IF (AMETHOD .GT. 0) CALL COMPUTE_SPATIAL()

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER SPATIAL ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN COMPUTE 4')
         ENDIF

#ifdef TIMINGS
         CALL MY_WTIME(TIME5)
#endif
         IF (SMETHOD .GT. 0 .AND. .NOT. (LSOURCESWAM .OR. LSOURCESWWIII)) THEN 
           CALL COMPUTE_SOURCES_EXP
         ELSE IF (SMETHOD .GT. 0 .AND. LSOURCESWAM) THEN
           DO IS = 1, MSC
             DO ID = 1, MDC
               FL3(:,ID,IS) =  AC2(:,IS,ID) * PI2 * SPSIG(IS)
             END DO
           END DO

           FL = FL3 
           THWOLD(:,1) = THWNEW
           U10NEW = MAX(TWO,SQRT(WINDXY(:,1)**2+WINDXY(:,2)**2)) * WINDFAC
           
           DO IP = 1, MNP
 
              THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))
!              WRITE(*,*) THWNEW(IP)
!              THWNEW(IP) = THWNEW(IP)*RADDEG
!              CALL DEG2NAUT (THWNEW(IP), DEG, .TRUE.) 
!              THWNEW(IP) = DEG
!              WRITE(*,*) THWNEW(IP)
!              STOP

              IF (ABS(IOBP(IP)) .GT. 0) CYCLE

         WRITE(111112,'(A10,I10)') 'BEFORE', IP
         WRITE(111112,'(A10,F20.10)') 'FL3', SUM(FL3(IP,:,:))
         WRITE(111112,'(A10,F20.10)') 'FL', SUM(FL(IP,:,:))
         WRITE(111112,'(A10,F20.10)') 'THWOLD', THWOLD(IP,1)
         WRITE(111112,'(A10,F20.10)') 'USOLD', USOLD(IP,1)
         WRITE(111112,'(A10,F20.10)') 'U10NEW', U10NEW(IP)
         WRITE(111112,'(A10,F20.10)') 'THWNEW', THWNEW(IP)
         WRITE(111112,'(A10,F20.10)') 'Z0OLD', Z0OLD(IP,1)
         WRITE(111112,'(A10,F20.10)') 'TAUW', TAUW(IP)
         WRITE(111112,'(A10,F20.10)') 'ROAIRO', ROAIRO(IP,1)
         WRITE(111112,'(A10,F20.10)') 'ZIDLOLD', ZIDLOLD(IP,1)
         WRITE(111112,'(A10,F20.10)') 'Z0NEW', Z0NEW(IP)
         WRITE(111112,'(A10,F20.10)') 'ROAIRN', ROAIRN(IP)
         WRITE(111112,'(A10,F20.10)') 'ZIDLNEW', ZIDLNEW(IP)
         WRITE(111112,'(A10,F20.10)') 'SL', SUM(SL(IP,:,:))
         WRITE(111112,'(A10,F20.10)') 'FCONST', SUM(FCONST(IP,:))

           CALL IMPLSCH (FL3(IP,:,:), FL(IP,:,:), IP, IP, 1, &
     &                   THWOLD(IP,1), USOLD(IP,1), &
     &                   TAUW(IP), Z0OLD(IP,1), &
     &                   ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                   U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                   Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                   SL(IP,:,:), FCONST(IP,:))

         WRITE(111112,'(A10,I10)') 'AFTER', IP 
         WRITE(111112,'(A10,F20.10)') 'FL3', SUM(FL3(IP,:,:))
         WRITE(111112,'(A10,F20.10)') 'FL', SUM(FL(IP,:,:))
         WRITE(111112,'(A10,F20.10)') 'THWOLD', THWOLD(IP,1)
         WRITE(111112,'(A10,F20.10)') 'USOLD', USOLD(IP,1)
         WRITE(111112,'(A10,F20.10)') 'U10NEW', U10NEW(IP)
         WRITE(111112,'(A10,F20.10)') 'THWNEW', THWNEW(IP)
         WRITE(111112,'(A10,F20.10)') 'Z0OLD', Z0OLD(IP,1)
         WRITE(111112,'(A10,F20.10)') 'TAUW', TAUW(IP)
         WRITE(111112,'(A10,F20.10)') 'ROAIRO', ROAIRO(IP,1)
         WRITE(111112,'(A10,F20.10)') 'ZIDLOLD', ZIDLOLD(IP,1)
         WRITE(111112,'(A10,F20.10)') 'Z0NEW', Z0NEW(IP)
         WRITE(111112,'(A10,F20.10)') 'ROAIRN', ROAIRN(IP)
         WRITE(111112,'(A10,F20.10)') 'ZIDLNEW', ZIDLNEW(IP)
         WRITE(111112,'(A10,F20.10)') 'SL', SUM(SL(IP,:,:))
         WRITE(111112,'(A10,F20.10)') 'FCONST', SUM(FCONST(IP,:))

           END DO

           DO IS = 1, MSC
             DO ID = 1, MDC
               AC2(:,IS,ID) =  FL3(:,ID,IS) / PI2 / SPSIG(IS)
             END DO
           END DO
         ELSE IF (SMETHOD .GT. 0 .AND. LSOURCESWWIII) THEN 
           !!!!
         ENDIF

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER SOURCES ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN COMPUTE 5')
         ENDIF

#ifdef TIMINGS
         CALL MY_WTIME(TIME6)
#endif

         IF (LMAXETOT .AND. SMETHOD .EQ. 0) CALL BREAK_LIMIT_ALL ! Miche for no source terms ... may cause oscilations ...

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER BREAK LIMIT ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN COMPUTE 6')
         ENDIF

#ifdef TIMINGS
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----SIMPLE SPLITTING SCHEME-----'
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'DIFFRACTION                      ', TIME2-TIME1
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SOURCES                          ', TIME6-TIME5
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS ADVEKTION            ', TIME5-TIME4
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS THETA SPACE          ', TIME4-TIME3
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS SIGMA SPACE          ', TIME3-TIME2
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU MICHE LIMITER                ', TIME6-TIME5
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS TOTAL TIME           ', TIME6-TIME1
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-------------TIMINGS-------------'
#endif
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE COMPUTE_SIMPLE_EXPLICIT'
         FLUSH(STAT%FHNDL)

        IF (.NOT. LDIFR) LCALC = .FALSE.

!        CALL PLOT_SHADED_CONTOUR_POLAR(SPSIG/PI2,SPDIR*RADDEG,MSC,MDC,AC2(137,:,:),10,MSC,MDC,'BEFORE ANY CALL')
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DOUBLE_STRANG_EXPLICIT()
        USE DATAPOOL
        IMPLICIT NONE

        REAL(rkind) :: TIME1, TIME2, TIME3, TIME4, TIME5
        REAL(rkind) :: TIME6, TIME7, TIME8, TIME9, TIME10
        REAL(rkind) :: TIME11, TIME12, TIME13, TIME14, TIME15, TIME16, TIME17

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'START COMPUTE'

        IF (.NOT. LSTEA .AND. .NOT. LQSTEA) THEN
          DT4A = 0.5*MAIN%DELT
          DT4S = 0.5*DT4A
          DT4D = ONETHIRD*MAIN%DELT
          DT4F = DT4D 
        ELSE IF (LQSTEA) THEN
          DT4A = 0.5*DT_ITER
          DT4S = DT4A * 0.25
          DT4D = ONETHIRD*DT_ITER
          DT4F = DT4D 
        END IF

#ifdef TIMINGS
        CALL MY_WTIME(TIME1)
#endif

        CALL COMPUTE_DIFFRACTION

        IF (FMETHOD .GT. 0) CALL COMPUTE_FREQUENCY
        IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION
! ---- 1st spectra 
        IF (AMETHOD .GT. 0) CALL COMPUTE_SPATIAL
        IF (SMETHOD .GT. 0) CALL COMPUTE_SOURCES_EXP
        IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION
        IF (FMETHOD .GT. 0) CALL COMPUTE_FREQUENCY
! ---- 2nd spectra 
        IF (AMETHOD .GT. 0) CALL COMPUTE_SPATIAL()
        IF (FMETHOD .GT. 0) CALL COMPUTE_FREQUENCY
        IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION()
! ---- 3rd spectra 
        IF (AMETHOD .GT. 0) CALL COMPUTE_SPATIAL()
        IF (SMETHOD .GT. 0) CALL COMPUTE_SOURCES_EXP()

#ifdef TIMINGS
        CALL MY_WTIME(TIME17)
#endif

        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '------DOUBLE STRANG SPLITTING----'
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'DIFFRACTION                      ', TIME2-TIME1 + TIME5-TIME4 + TIME8-TIME7 + TIME14+TIME13
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SOURCES                          ', TIME7-TIME6 + TIME13-TIME12
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS ADVEKTION            ', TIME6-TIME5 + TIME12-TIME11
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS THETA SPACE          ', TIME4-TIME3 + TIME9-TIME8 + TIME16-TIME15
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS SIGMA SPACE          ', TIME3-TIME2 + TIME10-TIME9 + TIME15-TIME14 
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS TOTAL TIME           ', TIME17-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-------------TIMINGS-------------'

        IF (.NOT. LDIFR) LCALC = .FALSE.
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_IMPLICIT
        USE DATAPOOL
        IMPLICIT NONE

        REAL(rkind), SAVE       :: TIME1, TIME2, TIME3, TIME4, TIME5, TIME6, TIME7, TIME8, TIME9, TIME10, TIME11
        INTEGER          :: IP, IT


        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE COMPUTE_IMPLICIT'
        FLUSH(STAT%FHNDL)

        IF (.NOT. LSTEA .AND. .NOT. LQSTEA) THEN
          DT4A = MAIN%DELT
          DT4S = DT4A
          DT4D = 0.5_rkind*DT4A
          DT4F = 0.5_rkind*DT4A
        ELSE IF (LQSTEA) THEN
          DT4A = DT_ITER
          DT4S = DT4A
          DT4D = 0.5_rkind*DT4A
          DT4F = 0.5_rkind*DT4A
        END IF

        AC1  = AC2

#ifdef TIMINGS
        CALL MY_WTIME(TIME1)
#endif

        CALL COMPUTE_DIFFRACTION

#ifdef TIMINGS
        CALL MY_WTIME(TIME2)
#endif
        IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION
        IF (FMETHOD .GT. 0) CALL COMPUTE_FREQUENCY
        IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION
        IF (FMETHOD .GT. 0) CALL COMPUTE_FREQUENCY
#ifdef TIMINGS
        CALL MY_WTIME(TIME3)
#endif
        IF (LMAXETOT) CALL BREAK_LIMIT_ALL ! Enforce Miche
#ifdef TIMINGS
        CALL MY_WTIME(TIME4)
#endif
        IF (SMETHOD .GT. 0) CALL COMPUTE_SOURCES_IMP
#ifdef TIMINGS
        CALL MY_WTIME(TIME5)
#endif
        IF (AMETHOD .GT. 0) CALL COMPUTE_SPATIAL
#ifdef TIMINGS
       CALL MY_WTIME(TIME6)
#endif
        IF (LLIMT .AND. SMETHOD .GT. 0) CALL ACTION_LIMITER
#ifdef TIMINGS
        CALL MY_WTIME(TIME7)
#endif
        IF (LMAXETOT) CALL BREAK_LIMIT_ALL ! Enforce Miche  
#ifdef TIMINGS
        CALL MY_WTIME(TIME8)
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----IMPLICIT SPLITTING SCHEME-----'
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'DIFFRACTION                      ', TIME2-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS ADVEKTION            ', TIME6-TIME5
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS SPECTRAL SPACE       ', TIME3-TIME2
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS SOURCES              ', TIME5-TIME4
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'ACTION LIMITER                   ', TIME7-TIME6
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'MICHE LIMITER                    ', TIME8-TIME7+TIME4-TIME3
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS TOTAL TIME           ', TIME8-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-------------TIMINGS-------------'
        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE COMPUTE_IMPLICIT'
        FLUSH(STAT%FHNDL)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_ITERATIVE_SPLITTING()
        USE DATAPOOL
        IMPLICIT NONE

        INTEGER          :: ITER, IS, ID, IP

#ifdef TIMINGS
        REAL(rkind), SAVE       :: TIME1, TIME2
#endif

#ifdef TIMINGS
         CALL MY_WTIME(TIME1)
#endif

         DT4A = MAIN%DELT 
         DT4S = DT4A
         DT4D = DT4A
         DT4F = DT4A

! Set DAC's to Zero ... 

#ifdef TIMINGS
         CALL MY_WTIME(TIME1)
#endif
         DAC_THE = 0.
         DAC_SIG = 0.
         DAC_SOU = 0.
         DAC_ADV = 0.

         AC1 = AC2

! 1st step ...

         IITERSPLIT = 0

         CALL COMPUTE_SPATIAL()
         CALL COMPUTE_FREQUENCY
         CALL COMPUTE_DIRECTION()
         CALL COMPUTE_SOURCES_EXP()

         WRITE(DBG%FHNDL,*) SUM(DAC_ADV), SUM(DAC_SOU), SUM(DAC_THE)

         AC2 = AC1

         IITERSPLIT = 1

! iteration ...

         CALL COMPUTE_SPATIAL()
         CALL COMPUTE_FREQUENCY
         CALL COMPUTE_DIRECTION()
         CALL COMPUTE_SOURCES_EXP()

         WRITE(DBG%FHNDL,*) SUM(DAC_ADV), SUM(DAC_SOU), SUM(DAC_THE)

         AC2 = AC1

! final step ... 

         CALL COMPUTE_SPATIAL()
         CALL COMPUTE_FREQUENCY
         CALL COMPUTE_DIRECTION()
         CALL COMPUTE_SOURCES_EXP()

#ifdef TIMINGS
         CALL MY_WTIME(TIME2)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DIFFRACTION
        USE DATAPOOL
        IF (LDIFR) THEN
          IF (IDIFFR == 1 ) THEN
            CALL DIFFRA_SIMPLE
          ELSE IF (IDIFFR == 2) THEN
            CALL DIFFRA_EXTENDED
          END IF
        END IF
      END SUBROUTINE COMPUTE_DIFFRACTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_SPATIAL()
        USE DATAPOOL
        IMPLICIT NONE

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE COMPUTE_SPATIAL'
        FLUSH(STAT%FHNDL)

        IF (DIMMODE == 1) THEN
          CALL COMPUTE_ADVECTION1D_QUICKEST_A
        ELSE IF (DIMMODE == 2) THEN
          IF (LVECTOR) THEN
            CALL FLUCT_3
          ELSE
            CALL FLUCT_1
! don't forget to uncomment FLUCT* in wwm_fluctsplit
!             IF(ICOMP == 0) THEN
!               CALL FLUCT_EXPLICIT()
!             ELSE IF(ICOMP == 1) THEN
!               CALL FLUCT_SEMIIMPLICIT()
!             ELSE IF(ICOMP == 2) THEN
!               CALL FLUCT_IMPLICIT()
!             ENDIF
          END IF
          IF ( ICOMP .GE. 1 .AND. (AMETHOD .EQ. 2 .OR. AMETHOD .EQ. 3 )) CALL RESCALE_SPECTRUM
        END IF

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE_SPATIAL'
        FLUSH(STAT%FHNDL)

      END SUBROUTINE COMPUTE_SPATIAL
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DIRECTION()
        USE DATAPOOL
        IMPLICIT NONE

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE_DIRECTION'
        FLUSH(STAT%FHNDL)
 
        IF (DMETHOD > 0) THEN
          IF (DMETHOD == 1) THEN
            CALL COMPUTE_DIRECTION_CNTG_A
          ELSE IF (DMETHOD == 2) THEN
            CALL COMPUTE_DIRECTION_QUICKEST_A
          ELSE IF (DMETHOD == 3) THEN
            CALL COMPUTE_DIRECTION_WENO_A
          ELSE IF (DMETHOD == 4) THEN
            CALL COMPUTE_DIRECTION_UPWIND_A
          ELSE IF (DMETHOD == 5) THEN
            CALL COMPUTE_DIRECTION_UPWIND_IMPLICIT
          END IF
        END IF

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE_DIRECTION'
        FLUSH(STAT%FHNDL)

        IF ( DMETHOD == 1) CALL RESCALE_SPECTRUM

      END SUBROUTINE COMPUTE_DIRECTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_FREQUENCY()
        USE DATAPOOL
        IMPLICIT NONE

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE_FREQUENCY'
        FLUSH(STAT%FHNDL)

        IF (FMETHOD == 1) CALL COMPUTE_FREQUENCY_QUICKEST_A
        IF (FMETHOD == 2) CALL COMPUTE_FREQUENCY_UPWIND_EXPLICIT
        IF (FMETHOD == 3) CALL COMPUTE_FREQUENCY_UPWIND_IMPLICIT

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE_FREQUENCY'
        FLUSH(STAT%FHNDL)

      END SUBROUTINE COMPUTE_FREQUENCY
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_SOURCES_EXP()
        USE DATAPOOL
        IMPLICIT NONE

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE_SOURCES_EXP'
        FLUSH(STAT%FHNDL)

        IF (ICOMP < 2 .AND. SMETHOD > 0) THEN
          CALL SOURCE_INT_EXP
        END IF

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE_SOURCES_EXP'
        FLUSH(STAT%FHNDL)

      END SUBROUTINE COMPUTE_SOURCES_EXP
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_SOURCES_IMP()
        USE DATAPOOL
        IMPLICIT NONE

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE_SOURCES_IMP'
        FLUSH(STAT%FHNDL)

        IF (ICOMP >= 2  .AND. SMETHOD > 0) THEN
          CALL SOURCE_INT_IMP()
        END IF

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE_SOURCES_IMP'
        FLUSH(STAT%FHNDL)

      END SUBROUTINE COMPUTE_SOURCES_IMP
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CFLSPEC()
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER              :: IP

         REAL(rkind)                 :: TMPCFLCAD(MNP), TMPCAD(MNP)
         REAL(rkind)                 :: TMPCFLCAS(MNP), TMPCAS(MNP)
         REAL(rkind)                 :: CAS(MSC,MDC), CAD(MSC,MDC)

         OPEN(310, FILE='cflcad.bin', FORM = 'UNFORMATTED', STATUS = 'UNKNOWN')
         OPEN(311, FILE='cflcas.bin', FORM = 'UNFORMATTED', STATUS = 'UNKNOWN')

         TMPCFLCAS = 0.
         TMPCFLCAD = 0.
         TMPCAS    = 0. 
         TMPCAD    = 0.

         DO IP = 1, MNP
           IF (DEP(IP) .GT. DMIN) THEN
             CALL PROPTHETA(IP,CAD)
             CALL PROPSIGMA(IP,CAS)
             TMPCAD(IP)    = MAXVAL(ABS(CAD))
! 0.5 since the directional and frequency intergration is split in two parts ....
             TMPCFLCAD(IP) = 0.5 * TMPCAD(IP)*MAIN%DELT/DDIR
             TMPCAS(IP)    = MAXVAL(ABS(CAS))
! absolute max. value ... lies on the secure side ... to do ...
             TMPCFLCAS(IP) = 0.5 * TMPCAS(IP)*MAIN%DELT/MINVAL(DS_INCR)
           ELSE
             CALL PROPTHETA(IP,CAD)
             CALL PROPSIGMA(IP,CAS)
             TMPCFLCAD(IP) = 0.
             TMPCAD(IP)    = 0.
             TMPCFLCAS(IP) = 0.
             TMPCAS(IP)    = 0.
           END IF
         END DO

         MAXCFLCAD = MAXVAL(TMPCAD)
         MAXCFLCAS = MAXVAL(TMPCAS)

         WRITE (310) SNGL(RTIME)
         WRITE (310) (SNGL(TMPCAD(IP)), SNGL(TMPCAD(IP)), SNGL(TMPCFLCAD(IP)), IP = 1, MNP)
         WRITE (311) SNGL(RTIME)
         WRITE (311) (SNGL(TMPCAS(IP)), SNGL(TMPCAS(IP)), SNGL(TMPCFLCAS(IP)), IP = 1, MNP)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef PETSC
      SUBROUTINE COMPUTE_FULL_IMPLICIT_PATANKAR
      USE DATAPOOL
      USE PETSC_BLOCK, ONLY : FREQ_SHIFT_IMPL, REFRACTION_IMPL, SOURCE_IMPL, EIMPS_PETSC_BLOCK
      IMPLICIT NONE
      IF (.NOT. LSTEA .AND. .NOT. LQSTEA) THEN
        DT4A = MAIN%DELT
        DT4S = DT4A
        DT4D = 0.5_rkind*DT4A
        DT4F = 0.5_rkind*DT4A 
      ELSE IF (LQSTEA) THEN
        DT4A = DT_ITER
        DT4S = DT4A
        DT4D = 0.5_rkind*DT4A
        DT4F = 0.5_rkind*DT4A
      END IF
      AC1  = AC2
      CALL COMPUTE_DIFFRACTION
      !
      ! Below is for debugging purpose only. 
      ! Only used if the refraction/freq_shift/source implicit
      ! are not selected
      !
      CALL SOURCE_INT_IMP()
      IF (SOURCE_IMPL .eqv. .FALSE.) THEN
        ! Do something clearly
      END IF
      IF (REFRACTION_IMPL .eqv. .FALSE.) THEN
        IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION()
      END IF
      IF (FREQ_SHIFT_IMPL .eqv. .FALSE.) THEN
        IF (FMETHOD .GT. 0) CALL COMPUTE_FREQUENCY()
      END IF
      !
      ! the Mother of all implicit computations.
      !
      CALL EIMPS_PETSC_BLOCK

      IF (LLIMT .AND. SMETHOD .GT. 0) CALL ACTION_LIMITER
      IF (LMAXETOT) CALL BREAK_LIMIT_ALL ! Enforce Miche  
      END SUBROUTINE
#endif
