#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INT_WAM(IP,DT,ACLOC)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP
         REAL(rkind), INTENT(IN) :: DT
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)

         INTEGER       :: IS, ID
         REAL(rkind)   :: NEWDAC
         REAL(rkind)   :: SSINL(MSC,MDC)
         REAL(rkind)   :: SSINE(MSC,MDC),DSSINE(MSC,MDC)
         REAL(rkind)   :: SSDS(MSC,MDC),DSSDS(MSC,MDC)
         REAL(rkind)   :: SSNL3(MSC,MDC),DSSNL3(MSC,MDC)
         REAL(rkind)   :: SSNL4(MSC,MDC),DSSNL4(MSC,MDC)
         REAL(rkind)   :: SSBR(MSC,MDC),DSSBR(MSC,MDC)
         REAL(rkind)   :: SSBRL(MSC,MDC)
         REAL(rkind)   :: SSBF(MSC,MDC),DSSBF(MSC,MDC)
         REAL(rkind)   :: IMATRA(MSC,MDC), IMATDA(MSC,MDC), ACOLD(MSC,MDC)

         IF (IOBP(IP) .NE. 0 .AND. .NOT. LSOUBOUND) THEN
           RETURN
         ELSE IF (LSOUBOUND .AND. IOBP(IP) .EQ. 2) THEN
           RETURN
         ENDIF

         ACOLD = ACLOC
         SSINL = ZERO
         SSINE = ZERO; DSSINE = ZERO
         SSDS  = ZERO; DSSDS  = ZERO
         SSNL3 = ZERO; DSSNL3 = ZERO
         SSNL4 = ZERO; DSSNL4 = ZERO
         SSBR  = ZERO; DSSBR  = ZERO   
         SSBF  = ZERO; DSSBF  = ZERO 
         IMATDA = ZERO; IMATRA = ZERO
         SSBRL = ZERO

         CALL DEEP_WATER(IP, ACLOC, IMATRA, IMATDA, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
         IF (ISHALLOW(IP)) CALL SHALLOW_WATER(IP, ACLOC, IMATRA, IMATDA, SSBR, DSSBR, SSBF, DSSBF, SSBRL, SSNL3, DSSNL3)
         DO IS = 1, MSC
           DO ID = 1, MDC
             NEWDAC = IMATRA(IS,ID) * DT / (ONE-DT*MIN(ZERO,IMATDA(IS,ID)))
!             write(*,*) NEWDAC / ( IMATRA(IS,ID) * DT )
             ACLOC(IS,ID) = MAX( ZERO, ACOLD(IS,ID) + NEWDAC )
           END DO
         END DO
         !CALL POST_INTEGRATION(IP,ACLOC)
         IF (LLIMT) CALL LIMITER(IP,ACOLD,ACLOC)

         WRITE(*,'(A20,6E20.10)') 'WAVE ACTION', SUM(ACLOC), MINVAL(ACLOC), MAXVAL(ACLOC)
         WRITE(*,'(A20,6E20.10)') 'LINEAR INPUT', SUM(SSINL), MINVAL(SSINL), MAXVAL(SSINL)
         WRITE(*,'(A20,6E20.10)') 'BREAKING LIMITER', SUM(SSBRL), MINVAL(SSBRL), MAXVAL(SSBRL)
         WRITE(*,'(A20,6E20.10)') 'EXP INPUT', SUM(SSINE), SUM(DSSINE), MINVAL(SSINE), MAXVAL(SSINE), MINVAL(DSSINE), MAXVAL(DSSINE)
         WRITE(*,'(A20,6E20.10)') 'WHITECAP', SUM(SSDS), SUM(DSSDS), MINVAL(SSDS), MAXVAL(SSDS), MINVAL(DSSDS), MAXVAL(DSSDS)
         WRITE(*,'(A20,6E20.10)') 'SNL4', SUM(SSNL4), SUM(DSSNL4), MINVAL(SSNL4), MAXVAL(SSNL4), MINVAL(DSSNL4), MAXVAL(DSSNL4)
         WRITE(*,'(A20,6E20.10)') 'BOTTOM FRICTION', SUM(SSBF), SUM(DSSBF), MINVAL(SSBF), MAXVAL(SSBF), MINVAL(DSSBF), MAXVAL(DSSBF)
         WRITE(*,'(A20,6E20.10)') 'BREAKING', SUM(SSBR), SUM(DSSBR), MINVAL(SSBR), MAXVAL(SSBR), MINVAL(DSSBR), MAXVAL(DSSBR)
         WRITE(*,'(A20,6E20.10)') 'TOTAL SOURCE TERMS', SUM(IMATRA), SUM(IMATDA), MINVAL(IMATRA), MAXVAL(IMATRA), MINVAL(IMATDA), MAXVAL(IMATDA)
         
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INT_PATANKAR(IP,ACLOC,IMATRA,IMATDA)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)      :: IP
         REAL(rkind), INTENT(IN)  :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(OUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         REAL(rkind)   :: SSINL(MSC,MDC)
         REAL(rkind)   :: SSINE(MSC,MDC),DSSINE(MSC,MDC)
         REAL(rkind)   :: SSDS(MSC,MDC),DSSDS(MSC,MDC)
         REAL(rkind)   :: SSNL3(MSC,MDC),DSSNL3(MSC,MDC)
         REAL(rkind)   :: SSNL4(MSC,MDC),DSSNL4(MSC,MDC)
         REAL(rkind)   :: SSBR(MSC,MDC),DSSBR(MSC,MDC)
         REAL(rkind)   :: SSBRL(MSC,MDC)
         REAL(rkind)   :: SSBF(MSC,MDC),DSSBF(MSC,MDC)
         REAL(rkind)   :: ACOLD(MSC,MDC)

         ACOLD = ACLOC
         SSINL = ZERO
         SSINE = ZERO; DSSINE = ZERO
         SSDS  = ZERO; DSSDS  = ZERO
         SSNL3 = ZERO; DSSNL3 = ZERO
         SSNL4 = ZERO; DSSNL4 = ZERO
         SSBR  = ZERO; DSSBR  = ZERO
         SSBF  = ZERO; DSSBF  = ZERO
         IMATDA = ZERO; IMATRA = ZERO
         SSBRL = ZERO

         IF (IOBP(IP) .NE. 0 .AND. .NOT. LSOUBOUND) THEN
           RETURN
         ELSE IF (LSOUBOUND .AND. IOBP(IP) .EQ. 2) THEN
           RETURN
         ENDIF

         CALL DEEP_WATER(IP, ACLOC, IMATRA, IMATDA, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
         !IF (ISHALLOW(IP)) CALL SHALLOW_WATER(IP, ACLOC, IMATRA, IMATDA, SSBR, DSSBR, SSBF, DSSBF, SSBRL, SSNL3, DSSNL3)

         WRITE(*,'(A20,6E20.10)') 'WAVE ACTION', SUM(ACLOC), MINVAL(ACLOC), MAXVAL(ACLOC)
         WRITE(*,'(A20,6E20.10)') 'LINEAR INPUT', SUM(SSINL), MINVAL(SSINL), MAXVAL(SSINL)
         WRITE(*,'(A20,6E20.10)') 'BREAKING LIMITER', SUM(SSBRL), MINVAL(SSBRL), MAXVAL(SSBRL)
         WRITE(*,'(A20,6E20.10)') 'EXP INPUT', SUM(SSINE), SUM(DSSINE), MINVAL(SSINE), MAXVAL(SSINE), MINVAL(DSSINE), MAXVAL(DSSINE)
         WRITE(*,'(A20,6E20.10)') 'WHITECAP', SUM(SSDS), SUM(DSSDS), MINVAL(SSDS), MAXVAL(SSDS), MINVAL(DSSDS), MAXVAL(DSSDS)
         WRITE(*,'(A20,6E20.10)') 'SNL4', SUM(SSNL4), SUM(DSSNL4), MINVAL(SSNL4), MAXVAL(SSNL4), MINVAL(DSSNL4), MAXVAL(DSSNL4)
         WRITE(*,'(A20,6E20.10)') 'BOTTOM FRICTION', SUM(SSBF), SUM(DSSBF), MINVAL(SSBF), MAXVAL(SSBF), MINVAL(DSSBF), MAXVAL(DSSBF)
         WRITE(*,'(A20,6E20.10)') 'BREAKING', SUM(SSBR), SUM(DSSBR), MINVAL(SSBR), MAXVAL(SSBR), MINVAL(DSSBR), MAXVAL(DSSBR)
         WRITE(*,'(A20,6E20.10)') 'TOTAL SOURCE TERMS', SUM(IMATRA), SUM(IMATDA), MINVAL(IMATRA), MAXVAL(IMATRA), MINVAL(IMATDA), MAXVAL(IMATDA)
   
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE POST_INTEGRATION(IP, ACLOC)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)       :: IP
         REAL(rkind),INTENT(INOUT) :: ACLOC(MSC,MDC)
         INTEGER                   :: IS, ID
         REAL(rkind)               :: VEC2RAD, WINDTH, FPM, WIND10
         REAL(rkind)               :: IMATRA(MSC,MDC)
         REAL(rkind)               :: SSINE(MSC,MDC),DSSINE(MSC,MDC), SSINL(MSC,MDC)
         REAL(rkind)               :: SSDS(MSC,MDC),DSSDS(MSC,MDC)

         IF (ISOURCE == 1) THEN
           CALL ST4_POST(IP, ACLOC, SSINE, DSSINE, SSDS, DSSDS, SSINL)
         ELSE IF (ISOURCE == 2) THEN
           CALL ECMWF_POST(IP, ACLOC, SSINE, DSSINE, SSDS, DSSDS, SSINL)
         ELSE IF (ISOURCE == 3) THEN
!2do write some post code for cycle3
         ENDIF
      END SUBROUTINE POST_INTEGRATION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SOURCES_EXPLICIT
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP
         REAL(rkind) :: ACLOC(MSC,MDC)

         DO IP = 1, MNP
           ACLOC = AC2(:,:,IP)
           IF (SMETHOD == 1) THEN
             CALL INT_WAM(IP,DT4S,ACLOC)
           ENDIF
           AC2(:,:,IP) = ACLOC 
         ENDDO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SOURCES_IMPLICIT
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER     :: IP
         REAL(rkind) :: ACLOC(MSC,MDC),IMATRA(MSC,MDC),IMATDA(MSC,MDC)

         DO IP = 1, MNP
           ACLOC = AC2(:,:,IP)
           IF (SMETHOD == 1) THEN
             CALL INT_PATANKAR(IP,ACLOC,IMATRA,IMATDA)
           ENDIF
           IMATRAA(:,:,IP) = IMATRA
           IMATDAA(:,:,IP) = IMATDA
         ENDDO

         WRITE(*,*) 'SOURCES_IMPLICIT', SUM(IMATRAA), SUM(IMATDAA)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
         SUBROUTINE LIMITER(IP,ACOLD,ACLOC)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IP
         INTEGER                    :: IS, ID
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(IN)    :: ACOLD(MSC,MDC)
         REAL(rkind)                :: NEWDAC, OLDAC, NEWAC, DELT, XIMP, DELFL(MSC)
         REAL(rkind)                :: MAXDAC, CONST, SND, UFR_LIM, DELT5, USFM
         REAL(rkind)                :: MAXDACOLD

         CONST  = PI2**2*3.0*1.0E-7*DT4S*SPSIG(MSC)
         SND    = PI2*5.6*1.0E-3
         DELT   = DT4S
         XIMP   = 1._rkind
         DELT5  = XIMP*DELT
         DELFL  = COFRM4*DELT
         MAXDAC = ZERO

         DO IS = 1, MSC
           IF (ISOURCE .EQ. 1) THEN
             IF (UFRIC(IP) .GT. SMALL) THEN
               USFM   = UFRIC(IP)*MAX(FMEANWS(IP),FMEAN(IP))
               MAXDAC = USFM*DELFL(IS)/PI2/SPSIG(IS)
             ELSE
               LIMFAK = 0.1
               MAXDAC = 0.0081*LIMFAK/(TWO*SPSIG(IS)*WK(IS,IP)**3*CG(IS,IP))
             END IF
           ELSE IF (ISOURCE .EQ. 2) THEN
             IF (USNEW(IP) .GT. SMALL) THEN
               USFM   = USNEW(IP)*MAX(FMEANWS(IP),FMEAN(IP))
               MAXDAC = USFM*DELFL(IS)/PI2/SPSIG(IS)
             ELSE
               LIMFAK = 0.1
               MAXDAC = 0.0081*LIMFAK/(TWO*SPSIG(IS)*WK(IS,IP)**3*CG(IS,IP))
             END IF
           ELSE IF (ISOURCE .EQ. 3) THEN
             LIMFAK = 0.1
             MAXDAC = 0.0081*LIMFAK/(TWO*SPSIG(IS)*WK(IS,IP)**3*CG(IS,IP))
           END IF
           DO ID = 1, MDC
             NEWAC  = ACLOC(IS,ID)
             OLDAC  = ACOLD(IS,ID)
             NEWDAC = NEWAC - OLDAC
             NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)), NEWDAC)
             ACLOC(IS,ID) = MAX( ZERO, OLDAC + NEWDAC )
           END DO
         END DO

         END SUBROUTINE LIMITER
!**********************************************************************
!*                                                                    *
!**********************************************************************
         SUBROUTINE BREAK_LIMIT(IP,ACLOC,SSBRL)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)  :: IP

         REAL(rkind), INTENT(INOUT)  :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(OUT)    :: SSBRL(MSC,MDC)

         REAL(rkind)                 :: EFTAIL, HS
         REAL(rkind)                 :: EMAX, RATIO, ETOT
         REAL(rkind)                 :: DINTSPEC

         ETOT   = 0.0
         EFTAIL = 1.0 / (PTAIL(1)-1.0)

         ETOT = DINTSPEC(IP,ACLOC)

         HS = 4.*SQRT(ETOT)

         !IF (LMONO_IN) HMAX(IP) = HMAX(IP) * SQRT(2.)

         EMAX = 1./16. * (HMAX(IP))**2

         IF (ETOT .GT. EMAX .AND. ETOT .GT. THR) THEN
           RATIO = EMAX/ETOT
           SSBRL = ACLOC - RATIO * ACLOC
           ACLOC = RATIO * ACLOC
         END IF
         END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
         SUBROUTINE RESCALE_SPECTRUM
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP, IS, ID
         REAL(rkind)    :: ETOTAL, EPOSIT
         REAL(rkind)    :: FACTOR
         LOGICAL :: ENEG

         DO IP = 1, MNP
            DO IS = 1, MSC
               ETOTAL = ZERO 
               EPOSIT = ZERO 
               ENEG   = .FALSE.
               DO ID = 1, MDC
                  ETOTAL = ETOTAL + AC2(IS,ID,IP)
                  IF (AC2(IS,ID,IP) > ZERO) THEN
                     EPOSIT = EPOSIT + AC2(IS,ID,IP)
                  ELSE
                     ENEG = .TRUE.
                  END IF
               END DO
               IF (ENEG) THEN
                  IF (EPOSIT .GT. VERYSMALL) THEN
                    FACTOR = ETOTAL/EPOSIT
                  ELSE 
                    FACTOR = 0.
                  END IF
                  DO ID = 1, MDC
                     IF (AC2(IS,ID,IP) < ZERO) AC2(IS,ID,IP) = ZERO 
                     IF (FACTOR >= ZERO)  AC2(IS,ID,IP) = AC2(IS,ID,IP)*FACTOR
                     AC2(IS,ID,IP) = MAX(zero,AC2(IS,ID,IP))
                  END DO
               END IF
            END DO
         END DO
         END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
         SUBROUTINE SETSHALLOW
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP

         DO IP = 1, MNP
           IF (WK(1,IP)*DEP(IP) .LT. PI) THEN
             ISHALLOW(IP) = 1
           ELSE
             ISHALLOW(IP) = 0
           END IF
         END DO
         END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
