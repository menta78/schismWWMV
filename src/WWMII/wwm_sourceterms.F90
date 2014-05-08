#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ACTION_LIMITER()
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER                 :: IP, IS, ID
         REAL(rkind)             :: NEWDAC, OLDAC, NEWAC, DELT, XIMP, DELFL(MSC)
         REAL(rkind)             :: MAXDAC, CONST, SND, UFR_LIM, DELT5, USFM
         REAL(rkind)             :: ACLOC(MSC,MDC), ACOLD(MSC,MDC), MAXDACOLD


         CONST = PI2**2*3.0*1.0E-7*DT4S*SPSIG(MSC)
         SND   = PI2*5.6*1.0E-3

         DELT = DT4S
         XIMP = 1._rkind
         DELT5 = XIMP*DELT
         DELFL= COFRM4*DELT
         MAXDAC = ZERO

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(IP,ACOLD,ACLOC, MAXDAC, UFR_LIM, NEWAC, OLDAC, NEWDAC)
         DO IP = 1, MNP
           IF (DEP(IP) .LT. DMIN .OR. IOBP(IP) .EQ. 2) CYCLE
           ACOLD = AC1(IP,:,:)
           ACLOC = AC2(IP,:,:)
           DO IS = 1, MSC
             IF (MELIM .EQ. 1) THEN
               MAXDAC = 0.0081*LIMFAK/(TWO*SPSIG(IS)*WK(IP,IS)**3*CG(IP,IS))
             ELSE IF (MELIM .EQ. 2) THEN
               UFR_LIM = MAX(UFRIC(IP),G9*SND/SPSIG(IS))
               MAXDAC  = LIMFAK*ABS((CONST*UFR_LIM)/(SPSIG(IS)**3*WK(IP,IS)))
             ELSE IF (MELIM .EQ. 3) THEN
               IF (USNEW(IP) .GT. SMALL) THEN
                 MAXDAC = COFRM4(IS)*USNEW(IP)*MAX(FMEANWS(IP),FMEAN(IP))/PI2/SPSIG(IS)*DT4A
               ELSE
                 MAXDAC = 0.0081*LIMFAK/(TWO*SPSIG(IS)*WK(IP,IS)**3*CG(IP,IS))
               ENDIF
             END IF
             DO ID = 1, MDC
               NEWAC  = ACLOC(IS,ID)
               OLDAC  = ACOLD(IS,ID)
               NEWDAC = NEWAC - OLDAC
!               IF (NEWDAC .GT. 0.) THEN
                 NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)), NEWDAC)
!               ELSE
!                 IF (QBLOCAL(IP) .LT. THR) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)), NEWDAC)
!               END IF
               AC2(IP,IS,ID) = MAX( zero, OLDAC + NEWDAC )
               !IF (MELIM .EQ. 3) FL3(IP,ID,IS) = AC2(IP,IS,ID) * PI2 * SPSIG(IS)
             END DO
           END DO
         END DO
         END SUBROUTINE
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

!AR: Fuckung lmono shit crap !
         !IF (LMONO_IN) HMAX(IP) = HMAX(IP) * SQRT(2.)

         EMAX = 1./16. * (HMAX(IP))**2

         IF (ETOT .GT. EMAX) THEN
           RATIO = EMAX/ETOT
           SSBRL = ACLOC - RATIO * ACLOC
           ACLOC = RATIO * ACLOC
         END IF
        END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE BREAK_LIMIT_ALL
      USE DATAPOOL
      IMPLICIT NONE

      INTEGER              :: IP
      REAL(rkind)          :: HS
      REAL(rkind)          :: EMAX, RATIO, ETOT
      REAL(rkind)          :: DINTSPEC
      REAL(rkind)          :: ACLOC(MSC, MDC)
!      Print *, 'Passing BREAK_LIMIT_ALL'

      DO IP = 1, MNP
        ACLOC = AC2(IP,:,:)
        IF (ISHALLOW(IP) .EQ. 0) CYCLE

        ETOT = DINTSPEC(IP,ACLOC)

        HS = 4.*SQRT(ETOT)

        EMAX = 1./16. * (HMAX(IP))**2
!        WRITE(300,*) 'IP=', IP, ' HMAX=', HMAX(IP), ' DEP=', DEP(IP)
!        WRITE(300,*) '   ', IP, ' EMAX=', EMAX, ' ETOT=', ETOT
!        WRITE(300,*) '   ', IP, ' HS=', HS, ' BRHD=', BRHD

        IF (ETOT .GT. EMAX) THEN
          WRITE(300,*) '   break XP=', XP(IP)
          RATIO = EMAX/ETOT
          AC2(IP,:,:) = RATIO * ACLOC(:,:)
          AC1(IP,:,:) = RATIO * ACLOC(:,:)
        END IF
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE RESCALE_SPECTRUM()
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP, IS, ID
         REAL(rkind)    :: ETOTAL, EPOSIT
         REAL(rkind)    :: FACTOR
         LOGICAL :: ENEG

         DO IP = 1, MNP
            DO IS = 1, MSC
               ETOTAL = 0.0
               EPOSIT = 0.0
               ENEG   = .FALSE.
               DO ID = 1, MDC
                  ETOTAL = ETOTAL + AC2(IP,IS,ID)
                  IF (AC2(IP,IS,ID) > 0.0) THEN
                     EPOSIT = EPOSIT + AC2(IP,IS,ID)
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
                     IF (AC2(IP,IS,ID) < 0.0) AC2(IP,IS,ID) = 0.0
                     IF (FACTOR >= 0.0)  AC2(IP,IS,ID) = AC2(IP,IS,ID)*FACTOR
                     AC2(IP,IS,ID) = MAX(zero,AC2(IP,IS,ID))
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
           IF (WK(IP,1)*DEP(IP) .LT. PI) THEN
             ISHALLOW(IP) = 1
           ELSE
             ISHALLOW(IP) = 0
           END IF
         END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
