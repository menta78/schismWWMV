#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ECMWF_PRE (IP, ACLOC, IMATRA, IMATDA, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IP

         REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: SSINE(MSC,MDC),DSSINE(MSC,MDC), SSINL(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: SSDS(MSC,MDC),DSSDS(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: SSNL4(MSC,MDC),DSSNL4(MSC,MDC)

         INTEGER      :: IS, ID, IMETHOD

         REAL(rkind)  :: VEC2RAD
         REAL(rkind)  :: ETOT,SME01,SME10,KME01,KMWAM,KMWAM2,HS,WIND10
         REAL(rkind)  :: ETAIL,EFTAIL,EMAX,LIMAC,NEWDAC,FPM,WINDTH,TEMP,GTEMP1
         REAL(rkind)  :: RATIO,LIMFAC,LIMDAC,GTEMP2,FLHAB,DELFL,USFM, NEWDACDT
         REAL(rkind)  :: MAXDAC, MAXDACDT, MAXDACDTDA, SC, SP, DNEWDACDTDA, JAC

         IMETHOD = 0 
 
         SSINE = ZERO; DSSINE = ZERO
         SSDS  = ZERO; DSSDS  = ZERO
         SSNL4 = ZERO; DSSNL4 = ZERO
         SSINL = ZERO

!AR: next step is to use here the local routines ...
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

         CALL PREINTRHS (FL3(IP,:,:), FL(IP,:,:), IP, IP, 1, &
     &                   THWOLD(IP,1), USOLD(IP,1), &
     &                   TAUW(IP), Z0OLD(IP,1), &
     &                   ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                   U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                   Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                   SL(IP,:,:), FCONST(IP,:), FMEANWS(IP), MIJ(IP), &
     &                   SSDS, DSSDS, SSINE, DSSINE, &
     &                   SSNL4, DSSNL4)

         DO ID = 1, MDC
           DO IS = 1, MSC 
             JAC = ONE/PI2/SPSIG(IS)
             IF (IMETHOD == 0) THEN
               IMATRA(IS,ID) = SL(IP,ID,IS)*JAC
               IMATDA(IS,ID) = ZERO 
             ELSE IF (IMETHOD == 1) THEN 
               IMATRA(IS,ID) = (SSINE(ID,IS)+SSDS(ID,IS)+SSNL4(ID,IS))*JAC
             ELSE IF (IMETHOD == 2) THEN
               IMATRA(IS,ID) = (SSINE(ID,IS)+SSNL4(ID,IS))*JAC
               IMATDA(IS,ID) = -TWO*DSSDS(ID,IS)
             ELSE IF (IMETHOD == 3) THEN
               IF (SSINE(ID,IS) .GT. ZERO) THEN
                 IMATRA(IS,ID) = SSINE(ID,IS)*JAC
               ELSE
                 IMATDA(IS,ID) = -TWO*DSSINE(ID,IS)
               ENDIF
             ELSE IF (IMETHOD == 4) THEN
               GTEMP1 = MAX((1.-DT4A*FL(IP,ID,IS)),1.)
               GTEMP2 = DT4A*SL(IP,ID,IS)/GTEMP1
               FLHAB  = ABS(GTEMP2)
               DELFL  = COFRM4(IS)*DT4A
               USFM   = USNEW(IP)*MAX(FMEANWS(IP),FMEAN(IP))
               TEMP   = USFM*DELFL
               FLHAB  = MIN(FLHAB,TEMP)
               IMATRA(IS,ID) = SIGN(FLHAB,GTEMP2)/DT4A*JAC
               LIMFAC        = MIN(ONE,ABS(SIGN(FLHAB,GTEMP2/DT4A))/MAX(SMALL,ABS(IMATRAA(IS,ID,IP))))
               IMATDA(IS,ID) = ZERO ! -FL(IP,ID,IS)
             ENDIF 
           ENDDO ! ID
         ENDDO ! IS

         IF (.NOT. LINID) THEN
           CALL SET_WIND( IP, WIND10, WINDTH )
           CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
           CALL SIN_LIN_CAV(IP,WINDTH,FPM,IMATRA,SSINL)
           IMATRA = IMATRA + SSINL
         ENDIF

         WRITE(*,'(A20,6E20.10)') 'LINEAR INPUT', SUM(SSINL), MINVAL(SSINL), MAXVAL(SSINL)
         WRITE(*,'(A20,6E20.10)') 'WAVE ACTION', SUM(ACLOC), MINVAL(ACLOC), MAXVAL(ACLOC)
         WRITE(*,'(A20,6E20.10)') 'EXP INPUT', SUM(SSINE), SUM(DSSINE), MINVAL(SSINE), MAXVAL(SSINE), MINVAL(DSSINE), MAXVAL(DSSINE)
         WRITE(*,'(A20,6E20.10)') 'WHITECAP', SUM(SSDS), SUM(DSSDS), MINVAL(SSDS), MAXVAL(SSDS), MINVAL(DSSDS), MAXVAL(DSSDS)
         WRITE(*,'(A20,6E20.10)') 'SNL4', SUM(SSNL4), SUM(DSSNL4), MINVAL(SSNL4), MAXVAL(SSNL4), MINVAL(DSSNL4), MAXVAL(DSSNL4)
         WRITE(*,'(A20,6E20.10)') 'TOTAL SOURCE TERMS', SUM(IMATRA), SUM(IMATDA), MINVAL(IMATRA), MAXVAL(IMATRA), MINVAL(IMATDA), MAXVAL(IMATDA)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ECMWF_POST
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER        :: IP, IS, ID
         REAL(rkind)    :: ACLOC(MSC,MDC), VEC2RAD
         REAL(rkind)    :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         THWOLD(IP,1) = THWNEW(IP)
         THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))
         U10NEW(IP) = MAX(TWO,SQRT(WINDXY(IP,1)**2+WINDXY(IP,2)**2)) * WINDFAC
         Z0NEW(IP) = Z0OLD(IP,1)
         DO IS = 1, MSC
           DO ID = 1, MDC
             FL3(IP,ID,IS) =  AC2(IS,ID,IP) * PI2 * SPSIG(IS)
             FL(IP,ID,IS)  =  IMATDAA(IS,ID,IP)
             SL(IP,ID,IS)  =  IMATRAA(IS,ID,IP) * PI2 * SPSIG(IS)
           END DO
         END DO
         CALL POSTINTRHS (FL3(IP,:,:), FL(IP,:,:), IP, IP, 1, &
     &                    THWOLD(IP,1), USOLD(IP,1), &
     &                    TAUW(IP), Z0OLD(IP,1), &
     &                    ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                    U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                    Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                    SL(IP,:,:), FCONST(IP,:), FMEANWS(IP), MIJ(IP))
         DO IS = 1, MSC
           DO ID = 1, MDC
             AC2(IS,ID,IP) = FL3(IP,ID,IS) / PI2 / SPSIG(IS)
           END DO
         END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
