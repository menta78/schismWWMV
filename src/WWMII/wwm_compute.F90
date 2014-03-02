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
         IF (FMETHOD .GT. 0) CALL COMPUTE_FREQUENCY
         IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER DIRECTION AND FREQUENCY -2- ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN COMPUTE 3') 
         ENDIF

#ifdef TIMINGS
         CALL MY_WTIME(TIME4)
#endif
         IF (AMETHOD .GT. 0) CALL COMPUTE_SPATIAL

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
           FL = FL3 
           THWOLD(:,1) = THWNEW
           U10NEW = MAX(TWO,SQRT(WINDXY(:,1)**2+WINDXY(:,2)**2)) * WINDFAC
           DO IP = 1, MNP
             DO IS = 1, MSC
               DO ID = 1, MDC
                 FL3(1,ID,IS) = AC2(IP,IS,ID) * PI2 * SPSIG(IS)
                 FL(1,ID,IS)  = FL3(1,ID,IS)
                 SL(1,ID,IS)  = FL(1,ID,IS)
               END DO
               Z0NEW(IP) = Z0OLD(IP,1) 
             END DO
             THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))
             IF (LOUTWAM .AND. IP == TESTNODE) THEN
               WRITE(111112,'(A10,I10)') 'AFTER', IP
               WRITE(111112,'(A10,F20.10)') 'FL3', SUM(FL3(1,:,:))
               WRITE(111112,'(A10,F20.10)') 'FL', SUM(FL(1,:,:))
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
               WRITE(111112,'(A10,F20.10)') 'SL', SUM(SL(1,:,:))
               WRITE(111112,'(A10,F20.10)') 'FCONST', SUM(FCONST(1,:))
             ENDIF
             IF (.FALSE.) THEN
               CALL IMPLSCH (FL3(1,:,:), FL(1,:,:), IP, IP, 1, &
     &                       THWOLD(IP,1), USOLD(IP,1), &
     &                       TAUW(IP), Z0OLD(IP,1), &
     &                       ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                       U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                       Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                       SL(1,:,:), FCONST(1,:))
             ELSE
               CALL PREINTRHS (FL3(1,:,:), FL(1,:,:), IP, IP, 1, &
     &                         THWOLD(IP,1), USOLD(IP,1), &
     &                         TAUW(IP), Z0OLD(IP,1), &
     &                         ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                         U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                         Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                         SL(1,:,:), FCONST(1,:), FMEANWS(IP), MIJ(IP))
             IF (LOUTWAM .AND. IP == TESTNODE) THEN
               WRITE(111112,'(A10,I10)') 'AFTER', IP
               WRITE(111112,'(A10,F20.10)') 'FL3', SUM(FL3(1,:,:))
               WRITE(111112,'(A10,F20.10)') 'FL', SUM(FL(1,:,:))
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
               WRITE(111112,'(A10,F20.10)') 'SL', SUM(SL(1,:,:))
               WRITE(111112,'(A10,F20.10)') 'FCONST', SUM(FCONST(1,:))
             ENDIF
             CALL INTSPECWAM (FL3(1,:,:), FL(1,:,:), IP, IP, 1, &
     &                        THWOLD(IP,1), USOLD(IP,1), &
     &                        TAUW(IP), Z0OLD(IP,1), &
     &                        ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                        U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                        Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                        SL(1,:,:), FCONST(1,:), FMEANWS(IP), MIJ(IP))
             CALL POSTINTRHS (FL3(1,:,:), FL(1,:,:), IP, IP, 1, &
     &                        THWOLD(IP,1), USOLD(IP,1), &
     &                        TAUW(IP), Z0OLD(IP,1), &
     &                        ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                        U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                        Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                        SL(1,:,:), FCONST(1,:), FMEANWS(IP), MIJ(IP))
             ENDIF
             DO IS = 1, MSC
               DO ID = 1, MDC
                 AC2(IP,IS,ID) =  FL3(1,ID,IS) / PI2 / SPSIG(IS)
               END DO
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
      SUBROUTINE COMPUTE_SEMI_IMPLICIT
        USE DATAPOOL
        IMPLICIT NONE

        REAL(rkind), SAVE       :: TIME1, TIME2, TIME3, TIME4, TIME5, TIME6, TIME7, TIME8, TIME9, TIME10, TIME11, GTEMP1, GTEMP2
        INTEGER          :: IP, IT, IS, ID


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

        WRITE(*,*) '1', SUM(AC2)

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
        IF (SMETHOD .GT. 0 .AND. LSOURCESWAM) THEN 
          CALL SOURCE_INT_IMP_WAM_PRE 
        ELSE IF (SMETHOD .GT. 0 .AND. .NOT. LSOURCESWAM) THEN
          CALL SOURCE_INT_IMP_WWM
        ENDIF 

        WRITE(*,*) '2', SUM(IMATRAA), SUM(IMATDAA), SUM(AC2)
#ifdef TIMINGS
        CALL MY_WTIME(TIME5)
#endif
        IF (AMETHOD .GT. 0) CALL COMPUTE_SPATIAL
#ifdef TIMINGS
        CALL MY_WTIME(TIME6)
#endif
        WRITE(*,*) '3', SUM(IMATRAA), SUM(IMATDAA), SUM(AC2)
        IF (LLIMT .AND. SMETHOD .GT. 0 .AND. .NOT. LSOURCESWAM) CALL ACTION_LIMITER
        WRITE(*,*) '4', SUM(IMATRAA), SUM(IMATDAA), SUM(AC2)
        !IF (SMETHOD .GT. 0 .AND. LSOURCESWAM) CALL SOURCE_INT_IMP_WAM_POST
        WRITE(*,*) '5', SUM(IMATRAA), SUM(IMATDAA), SUM(AC2)
        CALL ACTION_LIMITER
        WRITE(*,*) '6', SUM(IMATRAA), SUM(IMATDAA), SUM(AC2)
#ifdef TIMINGS
        CALL MY_WTIME(TIME7)
#endif
        IF (LMAXETOT) CALL BREAK_LIMIT_ALL ! Enforce Miche  
#ifdef TIMINGS
        WRITE(*,*) '7', SUM(IMATRAA), SUM(IMATDAA), SUM(AC2)
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
      SUBROUTINE COMPUTE_ITERATIVE_SPLITTING
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

#ifdef TIMINGS
         CALL MY_WTIME(TIME1)
#endif
         DAC_THE = 0.
         DAC_SIG = 0.
         DAC_SOU = 0.
         DAC_ADV = 0.

         AC1 = AC2

         IITERSPLIT = 0

         CALL COMPUTE_SPATIAL()
         CALL COMPUTE_FREQUENCY
         CALL COMPUTE_DIRECTION()
         CALL COMPUTE_SOURCES_EXP()

         WRITE(DBG%FHNDL,*) SUM(DAC_ADV), SUM(DAC_SOU), SUM(DAC_THE)

         AC2 = AC1

         IITERSPLIT = 1

         CALL COMPUTE_SPATIAL()
         CALL COMPUTE_FREQUENCY
         CALL COMPUTE_DIRECTION()
         CALL COMPUTE_SOURCES_EXP()

         WRITE(DBG%FHNDL,*) SUM(DAC_ADV), SUM(DAC_SOU), SUM(DAC_THE)

         AC2 = AC1

         CALL COMPUTE_SPATIAL
         CALL COMPUTE_FREQUENCY
         CALL COMPUTE_DIRECTION
         CALL COMPUTE_SOURCES_EXP

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
      SUBROUTINE COMPUTE_SPATIAL
        USE DATAPOOL
        IMPLICIT NONE

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE COMPUTE_SPATIAL'
        FLUSH(STAT%FHNDL)

        IF (DIMMODE == 1) THEN
          CALL COMPUTE_ADVECTION1D_QUICKEST_A
        ELSE IF (DIMMODE == 2) THEN
          IF(ICOMP == 0) THEN
            CALL FLUCT_EXPLICIT
          ELSE IF(ICOMP == 1) THEN
            CALL FLUCT_SEMIIMPLICIT
          ELSE IF(ICOMP == 2) THEN
            CALL FLUCT_IMPLICIT
          ENDIF
          IF ( ICOMP .GE. 1 .AND. (AMETHOD .EQ. 2 .OR. AMETHOD .EQ. 3 )) CALL RESCALE_SPECTRUM
        END IF

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE_SPATIAL'
        FLUSH(STAT%FHNDL)

      END SUBROUTINE COMPUTE_SPATIAL
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DIRECTION
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
      SUBROUTINE COMPUTE_FREQUENCY
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
      SUBROUTINE COMPUTE_SOURCES_EXP
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
      SUBROUTINE COMPUTE_IMPLICIT
        USE DATAPOOL
        USE PETSC_BLOCK, ONLY : FREQ_SHIFT_IMPL, REFRACTION_IMPL, SOURCE_IMPL, EIMPS_PETSC_BLOCK
        IMPLICIT NONE
#ifdef TIMINGS
        REAL(rkind)       :: TIME1, TIME2, TIME3, TIME4, TIME5
        REAL(rkind)       :: TIME6, TIME7, TIME8, TIME9, TIME10, TIME11, TIME12, TIME13
#endif

        IF (.NOT. LSTEA .AND. .NOT. LQSTEA) THEN
          DT4A = MAIN%DELT
          DT4S = DT4A
          DT4D = DT4A
          DT4F = DT4A 
        ELSE IF (LQSTEA) THEN
          DT4A = DT_ITER
          DT4S = DT4A
          DT4D = DT4A
          DT4F = DT4A
        END IF

        AC1  = AC2

#ifdef TIMINGS
        CALL MY_WTIME(TIME1)
#endif
        CALL COMPUTE_DIFFRACTION
#ifdef TIMINGS
        CALL MY_WTIME(TIME2)
#endif
        IF (LMAXETOT) CALL BREAK_LIMIT_ALL ! Enforce Miche
#ifdef TIMINGS
        CALL MY_WTIME(TIME3)
#endif
        IF (SMETHOD .GT. 0 .AND. LSOURCESWAM) THEN
          CALL SOURCE_INT_IMP_WAM_PRE
        ELSE IF (SMETHOD .GT. 0 .AND. .NOT. LSOURCESWAM) THEN
          CALL SOURCE_INT_IMP_WWM
        ENDIF
#ifdef TIMINGS
        CALL MY_WTIME(TIME4)
#endif
        CALL EIMPS_PETSC_BLOCK
#ifdef TIMINGS
        CALL MY_WTIME(TIME5)
#endif
        IF (LLIMT .AND. SMETHOD .GT. 0) CALL ACTION_LIMITER
        IF (SMETHOD .GT. 0 .AND. LSOURCESWAM) CALL SOURCE_INT_IMP_WAM_POST
#ifdef TIMINGS
        CALL MY_WTIME(TIME6)
#endif
        IF (LMAXETOT) CALL BREAK_LIMIT_ALL ! Enforce Miche  
        IF (LLIMT .AND. SMETHOD .GT. 0) CALL ACTION_LIMITER
        IF (LMAXETOT) CALL BREAK_LIMIT_ALL ! Enforce Miche  
#ifdef TIMINGS
        CALL MY_WTIME(TIME7)
#endif
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
