#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ECMWF_PRE (IP, ACLOC, IMATRA, IMATDA, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IP

         REAL(rkind), INTENT(INOUT)    :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: SSINE(MSC,MDC),DSSINE(MSC,MDC), SSINL(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: SSDS(MSC,MDC),DSSDS(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: SSNL4(MDC,MSC),DSSNL4(MDC,MSC)

         INTEGER      :: IS, ID

         REAL(rkind)  :: VEC2RAD
         REAL(rkind)  :: ETOT,SME01,SME10,KME01,KMWAM,KMWAM2,HS,WIND10
         REAL(rkind)  :: ETAIL,EFTAIL,EMAX,LIMAC,NEWDAC,FPM,WINDTH,TEMP,GTEMP1
         REAL(rkind)  :: RATIO,LIMFAC,LIMDAC,GTEMP2,FLHAB,DELFL,USFM, NEWDACDT
         REAL(rkind)  :: MAXDAC, MAXDACDT, MAXDACDTDA, SC, SP, DNEWDACDTDA, JAC

         SSINE = ZERO; DSSINE = ZERO
         SSDS  = ZERO; DSSDS  = ZERO
         SSNL4 = ZERO; DSSNL4 = ZERO
         SSINL = ZERO

         DO IS = 1, MSC
           DO ID = 1, MDC
             FL3(IP,ID,IS) = ACLOC(IS,ID) * PI2 * SPSIG(IS)
             FL(IP,ID,IS)  = FL3(IP,ID,IS)
             SL(IP,ID,IS)  = FL(IP,ID,IS)
           END DO
         END DO

         THWOLD(IP,1) = THWNEW(IP)
         U10NEW(IP) = MAX(TWO,SQRT(WINDXY(IP,1)**2+WINDXY(IP,2)**2))*WINDFAC
         Z0NEW(IP) = Z0OLD(IP,1)
         THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))

!         CALL PREINTRHS (FL3(IP,:,:), FL(IP,:,:), IP, IP, 1, &
!     &                   THWOLD(IP,1), USOLD(IP,1), &
!     &                   TAUW(IP), Z0OLD(IP,1), &
!     &                   ROAIRO(IP,1), ZIDLOLD(IP,1), &
!     &                   U10NEW(IP), THWNEW(IP), USNEW(IP), &
!     &                   Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
!     &                   SL(IP,:,:), FCONST(IP,:), FMEANWS(IP), MIJ(IP), &
!     &                   SSDS, DSSDS, SSINE, DSSINE, &
!     &                   SSNL4, DSSNL4)

         CALL WAM_PRE (IP, FL3(IP,:,:), FL(IP,:,:), SL(IP,:,:), SSDS, DSSDS, SSNL4, DSSNL4, SSINE, DSSINE)

         DO ID = 1, MDC
           DO IS = 1, MSC 
             JAC = ONE/PI2/SPSIG(IS)
             IMATRA(IS,ID) = SL(IP,ID,IS)*JAC
             IMATDA(IS,ID) = FL(IP,ID,IS)/PI2
           ENDDO ! ID
         ENDDO ! IS

         IF (.NOT. LINID) THEN
           CALL SET_WIND( IP, WIND10, WINDTH )
           CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
           CALL SIN_LIN_CAV(IP,WINDTH,FPM,IMATRA,SSINL)
           IMATRA = IMATRA + SSINL
         ENDIF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ECMWF_POST(IP,ACLOC,SSINE,DSSINE,SSDS,DSSDS,SSINL)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)           :: IP
         REAL(rkind), INTENT(INOUT)    :: ACLOC(MSC,MDC)

         INTEGER                       :: IS, ID
         REAL(rkind)                   :: VEC2RAD, WINDTH, FPM, WIND10
         REAL(rkind)                   :: IMATRA(MSC,MDC)

         REAL(rkind), INTENT(OUT)   :: SSINE(MSC,MDC),DSSINE(MSC,MDC), SSINL(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: SSDS(MSC,MDC),DSSDS(MSC,MDC)

         THWOLD(IP,1) = THWNEW(IP)
         THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))
         U10NEW(IP) = MAX(TWO,SQRT(WINDXY(IP,1)**2+WINDXY(IP,2)**2)) * WINDFAC
         Z0NEW(IP) = Z0OLD(IP,1)
         DO IS = 1, MSC
           DO ID = 1, MDC
             FL3(IP,ID,IS) =  ACLOC(IS,ID) * PI2 * SPSIG(IS)
           END DO
         END DO
         CALL WAM_POST (IP, FL3(IP,:,:))
         DO IS = 1, MSC
           DO ID = 1, MDC
             ACLOC(IS,ID) = FL3(IP,ID,IS) / PI2 / SPSIG(IS)
           END DO
         END DO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
