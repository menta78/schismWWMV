#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SHALLOW_WATER (IP, ACLOC, IMATRA, IMATDA, SSBR, DSSBR, SSBF, DSSBF, SSBRL, SSNL3, DSSNL3)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)      :: IP
         REAL(rkind), INTENT(IN)  :: ACLOC(MSC,MDC)

         REAL(rkind), INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind), INTENT(OUT) :: SSNL3(MSC,MDC),DSSNL3(MSC,MDC)
         REAL(rkind), INTENT(OUT) :: SSBR(MSC,MDC),DSSBR(MSC,MDC)
         REAL(rkind), INTENT(OUT) :: SSBRL(MSC,MDC)
         REAL(rkind), INTENT(OUT) :: SSBF(MSC,MDC),DSSBF(MSC,MDC)

         INTEGER      :: IS, ID
         REAL(rkind)  :: VEC2RAD
         REAL(rkind)  :: ETOT,SME01,SME10,KME01,KMWAM,KMWAM2,HS
         REAL(rkind)  :: NEWDAC,FPM,WINDTH,TEMP
         REAL(rkind)  :: RATIO,LIMDAC,NEWDACDT
         REAL(rkind)  :: MAXDAC, SC, SP

         SSNL3 = ZERO; DSSNL3 = ZERO
         SSBR  = ZERO; DSSBR  = ZERO
         SSBRL = ZERO
         SSBF  = ZERO; DSSBF  = ZERO

         CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)

         IF (MESTR .GT. 0) CALL TRIAD_ELDEBERKY(IP, HS, SME01, ACLOC, SSNL3, DSSNL3)
         IF (MESBR .GT. 0) CALL SDS_SWB(IP, SME01, KMWAM, ETOT, HS, ACLOC, SSBR, DSSBR)
         IF (MESBF .GT. 0) CALL SDS_BOTF(IP,ACLOC,SSBF,DSSBF)

         IF (LSOURCESLIM) THEN
           DO IS = 1, MSC
             MAXDAC = 0.00081_rkind/(TWO*SPSIG(IS)*WK(IS,IP)**3*CG(IS,IP))
             DO ID = 1, MDC
               NEWDAC = SSNL3(IS,ID)*DT4A/(1.0-DT4A*MIN(ZERO,DSSNL3(IS,ID)))
               LIMDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
               SSNL3(IS,ID) = LIMDAC/DT4A
             ENDDO
           ENDDO
         ENDIF

         IF (ICOMP .GE. 2) THEN
           IF (optionCall .EQ. 1) THEN
             IMATDA = IMATDA + DSSBR + DSSNL3 + DSSBF
             IMATRA = IMATRA + SSBR  + SSNL3  + SSBF
           ELSE IF (optionCall .EQ. 2) THEN 
             IMATDA = IMATDA + DSSBR + DSSNL3 + DSSBF
             IMATRA = IMATRA + SSBR  + SSNL3  + SSBF
           ENDIF
         ELSE 
           IMATDA = IMATDA + DSSBR + DSSNL3 + DSSBF
           IMATRA = IMATRA + SSBR  + SSNL3  + SSBF
         ENDIF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
