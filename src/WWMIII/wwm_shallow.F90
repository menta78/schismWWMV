#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SHALLOW (IP, ACLOC, IMATRA, IMATDA)
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

         CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)

         IF (MESTR .GT. 0) CALL TRIAD_ELDEBERKY(IP, HS, SME01, ACLOC, IMATRA, IMATDA, SSNL3, DSSNL3)
         IF (MESBR .GT. 0) CALL SDS_SWB(IP, SME01, KMWAM, ETOT, HS, ACLOC, IMATRA, IMATDA, SSBR, DSSBR)
         IF (MESBF .GT. 0) CALL SDS_BOTF(IP,ACLOC,IMATRA,IMATDA,SSBF,DSSBF)

         IMATDA = IMATDA + DSSBR + DSSNL3 + DSSBF
         IMATRA = IMATRA + SSBR + SSNL3 

         IF (LMAXETOT) THEN
           NEWAC = ACLOC + IMATRA*DT4A/MAX((ONE-DT4A*IMATDA),ONE)
           CALL MEAN_WAVE_PARAMETER(IP,NEWAC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)
           EMAX = 1._rkind/16._rkind * (HMAX(IP))**2 ! HMAX is defined in the breaking routine or has some default value
           IF (ETOT .GT. EMAX) THEN
             RATIO  = EMAX/ETOT
             !SSBRL  = ACLOC*(RATIO-ONE)/DT4A
             DSSBRL = (RATIO-ONE)/DT4A ! Always negative ...
           END IF
           IMATDA = IMATDA + DSSBRL
         ENDIF

         WRITE(*,'(A20,6E20.10)') 'WAVE ACTION', SUM(ACLOC), MINVAL(ACLOC), MAXVAL(ACLOC)
         WRITE(*,'(A20,6E20.10)') 'BOTTOM FRICTION', SUM(SSBF), SUM(DSSBF), MINVAL(SSBF), MAXVAL(SSBF), MINVAL(DSSBF), MAXVAL(DSSBF)
         WRITE(*,'(A20,6E20.10)') 'BREAKING', SUM(SSBR), SUM(DSSBR), MINVAL(SSBR), MAXVAL(SSBR), MINVAL(DSSBR), MAXVAL(DSSBR)
         WRITE(*,'(A20,6E20.10)') 'BREAKING LIMITER', SUM(SSBRL), SUM(DSSBRL), MINVAL(SSBRL), MAXVAL(SSBRL), MINVAL(DSSBRL), MAXVAL(DSSBRL)
         !WRITE(*,'(A20,6E20.10)') 'LIMITER',  SUM(SSLIM), SUM(DSSLIM), MINVAL(SSLIM), MAXVAL(SSLIM), MINVAL(DSSLIM), MAXVAL(DSSLIM)
         WRITE(*,'(A20,6E20.10)') 'TOTAL SOURCE TERMS', SUM(IMATRA), SUM(IMATDA), MINVAL(IMATRA), MAXVAL(IMATRA), MINVAL(IMATDA), MAXVAL(IMATDA)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
