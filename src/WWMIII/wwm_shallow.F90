#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SHALLOW_WATER (IP, ACLOC, IMATRA, IMATDA)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IP

         REAL(rkind), INTENT(OUT)   :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)

         INTEGER      :: IS, ID, IMETHOD
         REAL(rkind)  :: VEC2RAD
         REAL(rkind)  :: NEWAC(MSC,MDC)
         REAL(rkind)  :: SSNL3(MSC,MDC),DSSNL3(MSC,MDC)
         REAL(rkind)  :: SSBR(MSC,MDC),DSSBR(MSC,MDC)
         REAL(rkind)  :: SSBRL(MSC,MDC),DSSBRL(MSC,MDC)
         REAL(rkind)  :: SSBF(MSC,MDC),DSSBF(MSC,MDC), SSINL(MSC,MDC)
         REAL(rkind)  :: ETOT,SME01,SME10,KME01,KMWAM,KMWAM2,HS,WIND10
         REAL(rkind)  :: ETAIL,EFTAIL,EMAX,LIMAC,NEWDAC,FPM,WINDTH,TEMP,GTEMP1
         REAL(rkind)  :: RATIO,LIMFAC,LIMDAC,GTEMP2,FLHAB,DELFL,USFM, NEWDACDT
         REAL(rkind)  :: MAXDAC, MAXDACDT, MAXDACDTDA, SC, SP, DNEWDACDTDA, JAC

         REAL(rkind),DIMENSION(MDC,MSC) :: SSDS,DSSDS,SSNL4,DSSNL4,SSINE,DSSINE

         SSNL3 = ZERO; DSSNL3 = ZERO
         SSBR  = ZERO; DSSBR  = ZERO
         SSBRL  = ZERO; DSSBRL  = ZERO
         SSBF  = ZERO; DSSBF  = ZERO

         IF (LMAXETOT) CALL BREAK_LIMIT(IP,ACLOC,SSBRL)

         CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)

         IF (MESTR .GT. 0) CALL TRIAD_ELDEBERKY(IP, HS, SME01, ACLOC, IMATRA, IMATDA, SSNL3, DSSNL3)
         IF (MESBR .GT. 0) CALL SDS_SWB(IP, SME01, KMWAM, ETOT, HS, ACLOC, IMATRA, IMATDA, SSBR, DSSBR)
         IF (MESBF .GT. 0) CALL SDS_BOTF(IP,ACLOC,IMATRA,IMATDA,SSBF,DSSBF)

         IF (LLIMT) THEN
           DO IS = 1, MSC
             MAXDAC = 0.00081_rkind/(TWO*SPSIG(IS)*WK(IS,IP)**3*CG(IS,IP))
             DO ID = 1, MDC
               NEWDAC = SSNL3(IS,ID)*DT4A/(1.0-DT4A*MIN(ZERO,DSSNL3(IS,ID)))
               LIMDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
               SSNL3(IS,ID) = LIMDAC/DT4A
             ENDDO
           ENDDO
         ENDIF

         IMATDA = IMATDA + DSSBR + DSSNL3 + DSSBF
         IMATRA = IMATRA + SSBR  + SSNL3 

         WRITE(*,'(A20,6E20.10)') 'WAVE ACTION', SUM(ACLOC), MINVAL(ACLOC), MAXVAL(ACLOC)
         WRITE(*,'(A20,6E20.10)') 'BOTTOM FRICTION', SUM(SSBF), SUM(DSSBF), MINVAL(SSBF), MAXVAL(SSBF), MINVAL(DSSBF), MAXVAL(DSSBF)
         WRITE(*,'(A20,6E20.10)') 'BREAKING', SUM(SSBR), SUM(DSSBR), MINVAL(SSBR), MAXVAL(SSBR), MINVAL(DSSBR), MAXVAL(DSSBR)
         WRITE(*,'(A20,6E20.10)') 'BREAKING LIMITER', SUM(SSBRL), SUM(DSSBRL), MINVAL(SSBRL), MAXVAL(SSBRL), MINVAL(DSSBRL), MAXVAL(DSSBRL)
         WRITE(*,'(A20,6E20.10)') 'TOTAL SOURCE TERMS', SUM(IMATRA), SUM(IMATDA), MINVAL(IMATRA), MAXVAL(IMATRA), MINVAL(IMATDA), MAXVAL(IMATDA)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
