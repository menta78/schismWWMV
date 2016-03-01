#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SCHALLOW (IP, ACLOC, IMATRA, IMATDA)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IP

         REAL(rkind), INTENT(OUT)   :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)

         INTEGER      :: IS, ID, IMETHOD
         REAL(rkind)  :: VEC2RAD
         REAL(rkind)  :: SSNL3(MSC,MDC),DSSNL3(MSC,MDC)
         REAL(rkind)  :: SSBR(MSC,MDC),DSSBR(MSC,MDC)
         REAL(rkind)  :: SSBF(MSC,MDC),DSSBF(MSC,MDC), SSINL(MSC,MDC)
         REAL(rkind)  :: ETOT,SME01,SME10,KME01,KMWAM,KMWAM2,HS,WIND10
         REAL(rkind)  :: ETAIL,EFTAIL,EMAX,LIMAC,NEWDAC,FPM,WINDTH,TEMP,GTEMP1
         REAL(rkind)  :: RATIO,LIMFAC,LIMDAC,GTEMP2,FLHAB,DELFL,USFM, NEWDACDT
         REAL(rkind)  :: MAXDAC, MAXDACDT, MAXDACDTDA, SC, SP, DNEWDACDTDA, JAC

         REAL(rkind),DIMENSION(MDC,MSC) :: SSDS,DSSDS,SSNL4,DSSNL4,SSINE,DSSINE

         IMETHOD = 4 

         SSNL3 = ZERO; DSSNL3 = ZERO
         SSBR  = ZERO; DSSBR  = ZERO
         SSBF  = ZERO; DSSBF  = ZERO
         SSINL = ZERO

         IF (IOBP(IP) .EQ. 0) THEN
           IF (ISHALLOW(IP) .EQ. 1) THEN
             CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)
             IF (MESTR .GT. 0) THEN
               CALL TRIAD_ELDEBERKY(IP, HS, SME01, ACLOC, IMATRA, IMATDA, SSNL3, DSSNL3)
               DO IS = 1, MSC
                 DO ID = 1, MDC
                   NEWDAC = SSNL3(IS,ID)*DT4A/MAX((1.-DT4A*DSSNL3(IS,ID)),1.)
                   MAXDAC = 0.0081*LIMFAK/(TWO*SPSIG(IS)*WK(IS,IP)**3*CG(IS,IP))*100
                   LIMFAC = ONE/MAX(ONE,NEWDAC/MAXDAC)
                   SC = SIGN(MIN(ABS(NEWDAC),MAXDAC),NEWDAC)/DT4A
                   !SSNL3(IS,ID)  = SC
                   !DSSNL3(IS,ID) = DSSNL3(IS,ID)*LIMFAC
                   !IF (ABS(SC) .GT. THR) WRITE(*,'(2I10,5F20.8)') IS, ID, NEWDAC, MAXDAC, DSSNL3(IS,ID), LIMFAC
                 END DO
               END DO
             ENDIF ! MESTR
             IF (MESBR .GT. 0) CALL SDS_SWB(IP, SME01, KMWAM, ETOT, HS, ACLOC, IMATRA, IMATDA, SSBR, DSSBR)
             IF (MESBF .GT. 0) CALL SDS_BOTF(IP,ACLOC,IMATRA,IMATDA,SSBF,DSSBF)
             IMATDA = IMATDA + DSSBR  + DSSNL3 + DSSBF
             IMATRA = IMATRA + SSBR + SSNL3 
           ENDIF ! ISHALLOW(IP) .EQ. 1
         ELSE ! IOBP(IP) .NE. 0
           IF (LSOUBOUND) THEN ! Source terms on boundary ...
             IF (IOBP(IP) .NE. 2) THEN
               IF (ISHALLOW(IP) .EQ. 1) THEN
                 CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)
                 IF (MESTR .GT. 0) THEN
                   CALL TRIAD_ELDEBERKY(IP, HS, SME01, ACLOC, IMATRA, IMATDA, SSNL3, DSSNL3)
                   DO IS = 1, MSC
                     DO ID = 1, MDC
                       NEWDAC = SSNL3(IS,ID)*DT4A/MAX((1.-DT4A*DSSNL3(IS,ID)),1.)
                       MAXDAC = 0.0081*LIMFAK/(TWO*SPSIG(IS)*WK(IS,IP)**3*CG(IS,IP))*100
                       LIMFAC = ONE/MAX(ONE,NEWDAC/MAXDAC)
                       SC = SIGN(MIN(ABS(NEWDAC),MAXDAC),NEWDAC)/DT4A
                       !SSNL3(IS,ID)  = SC
                       !DSSNL3(IS,ID) = DSSNL3(IS,ID)*LIMFAC
                       !IF (ABS(SC) .GT. THR) WRITE(*,'(2I10,5F20.8)') IS, ID, NEWDAC, MAXDAC, DSSNL3(IS,ID), LIMFAC
                     END DO
                   END DO
                 ENDIF ! MESTR
                 IF (MESBR .GT. 0) CALL SDS_SWB(IP, SME01, KMWAM, ETOT, HS, ACLOC, IMATRA, IMATDA, SSBR, DSSBR)
                 IF (MESBF .GT. 0) CALL SDS_BOTF(IP,ACLOC,IMATRA,IMATDA,SSBF,DSSBF)
                 IMATDA = IMATDA + DSSBR  + DSSNL3 + DSSBF
                 IMATRA = IMATRA + SSBR + SSNL3
               ENDIF ! ISHALLOW(IP) .EQ. 1
             ENDIF
           ENDIF
         ENDIF

         IF (IP == TESTNODE) THEN
           WRITE(*,'(A20,6E20.10)') 'WAVE ACTION', SUM(ACLOC), MINVAL(ACLOC), MAXVAL(ACLOC)
           WRITE(*,'(A20,6E20.10)') 'BOTTOM FRICTION', SUM(SSBF), SUM(DSSBF), MINVAL(SSBF), MAXVAL(SSBF), MINVAL(DSSBF), MAXVAL(DSSBF)
           WRITE(*,'(A20,6E20.10)') 'BREAKING', SUM(SSBR), SUM(DSSBR), MINVAL(SSBR), MAXVAL(SSBR), MINVAL(DSSBR), MAXVAL(DSSBR)
!           WRITE(*,'(A20,6E20.10)') 'BREAKING LIMITER', SUM(SSBRL), SUM(DSSBRL), MINVAL(SSBRL), MAXVAL(SSBRL), MINVAL(DSSBRL), MAXVAL(DSSBRL)
!           WRITE(*,'(A20,6E20.10)') 'LIMITER',  SUM(SSLIM), SUM(DSSLIM), MINVAL(SSLIM), MAXVAL(SSLIM), MINVAL(DSSLIM), MAXVAL(DSSLIM)
           WRITE(*,'(A20,6E20.10)') 'TOTAL SOURCE TERMS', SUM(IMATRA), SUM(IMATDA), MINVAL(IMATRA), MAXVAL(IMATRA), MINVAL(IMATDA), MAXVAL(IMATDA)
!           PAUSE
         ENDIF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
