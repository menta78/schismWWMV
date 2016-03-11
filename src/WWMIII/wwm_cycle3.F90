#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CYCLE3_PRE (IP, ACLOC, IMATRA, IMATDA, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IP
         REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)

         REAL(rkind), INTENT(OUT)   :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: SSINL(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: SSINE(MSC,MDC),DSSINE(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: SSDS(MSC,MDC),DSSDS(MSC,MDC)

         INTEGER                    :: IS, ID

         REAL(rkind)                :: NEWAC(MSC,MDC)
         REAL(rkind)                :: SSLIM(MSC,MDC), DSSLIM(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: SSNL4(MSC,MDC),DSSNL4(MSC,MDC)
         REAL(rkind)                :: ETOT,SME01,SME10,KME01,KMWAM
         REAL(rkind)                :: KMWAM2,HS,WIND10
         REAL(rkind)                :: EFTAIL,EMAX,NEWDAC,MAXDAC,FPM,WINDTH
         REAL(rkind)                :: RATIO,LIMDAC

         NEWAC = ZERO
         SSINL = ZERO
         SSINE = ZERO; DSSINE = ZERO
         SSNL4 = ZERO; DSSNL4 = ZERO
         SSDS   = ZERO; DSSDS = ZERO
         SSLIM  = ZERO; DSSLIM = ZERO
         IMATRA = ZERO; IMATDA = ZERO

         CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2) 

         CALL SET_WIND( IP, WIND10, WINDTH )
         CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
         IF (WIND10 .GT. THR .AND. ETOT .LT. THR) CALL SIN_LIN( IP, WINDTH, FPM, SSINL)
         IF (MESIN .GT. 0) CALL SIN_EXP( IP, WINDTH, ACLOC, SSINE, DSSINE )
         IF (MESDS .GT. 0) CALL SDS_CYCLE3_NEW ( IP, KMWAM, SME10, ETOT, ACLOC, SSDS, DSSDS )
         IF (MESNL .GT. 0) CALL SNL41(IP, KMWAM, ACLOC, IMATRA, IMATDA, SSNL4, DSSNL4)

         IMATRA = SSINL + SSDS + SSINE + SSNL4
         IMATDA = DSSDS + DSSNL4

         IF (LSOURCESLIM) THEN
           DO IS = 1, MSC
             MAXDAC = 0.00081_rkind/(TWO*SPSIG(IS)*WK(IS,IP)**3*CG(IS,IP))
             DO ID = 1, MDC
               NEWDAC = IMATRA(IS,ID)*DT4A/(1.0-DT4A*MIN(ZERO,IMATDA(IS,ID)))
               LIMDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
               IMATRA(IS,ID) = LIMDAC/DT4A
               !SSLIM(IS,ID)  = SIGN(ABS(NEWDAC-LIMDAC)/DT4A,NEWDAC)
               !DSSLIM(IS,ID) = SIGN(ABS(IMATDA(IS,ID)-ABS(LIMFAC*IMATDA(IS,ID))),NEWDAC)
             ENDDO
           ENDDO
         ENDIF

         WRITE(*,*) 'wind, th, fpm and sum sources'
         WRITE(*,*) WIND10, WINDTH, FPM, SUM(IMATRA), SUM(IMATDA)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDS_CYCLE3_NEW( IP, KMESPC, SMESPC, ETOT, ACLOC, SSDS, DSSDS )
!
!     Cycle 3 dissipation 
!
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)           :: IP
         REAL(rkind)   , INTENT(IN)    :: KMESPC, SMESPC, ETOT
         REAL(rkind)   , INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind)   , INTENT(OUT)   :: SSDS(MSC,MDC), DSSDS(MSC,MDC)

         INTEGER       :: IS

         REAL(rkind)   :: CDS, ALPHA_PM, FAC
         REAL(rkind)   :: STP_OV, STP_PM, N2
!
         ALPHA_PM  =  3.02E-3
         CDS       =  2.36E-5

         STP_OV = KMESPC * SQRT(ETOT)
         STP_PM = SQRT(ALPHA_PM)
         N2     = 4
         FAC    = CDS * (STP_OV / STP_PM)**N2

         DO IS = 1, MSC
           DSSDS(IS,:) = - FAC * SMESPC * (WK(IS,IP)/KMESPC)
           SSDS(IS,:)  =   DSSDS(IS,:) * ACLOC(IS,:)
         END DO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SIN_LIN( IP, WINDTH, FPM, SSINL )
!
!     Linear growth term according to Cavaleri & Melanotte Rizolli ...
!
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)   :: IP
         REAL(rkind)   , INTENT(OUT)  :: SSINL(MSC,MDC)
         REAL(rkind)   , INTENT(IN)   :: WINDTH
         REAL(rkind)   , INTENT(IN)   :: FPM

         INTEGER                      :: IS, ID
         REAL(rkind)                  :: AUX, AUX1, AUX2, AUXH
         REAL(rkind)                  :: SWINA

         AUX = 0.0015_rkind / ( G9*G9*PI2 )

         DO IS = 1, MSC
           AUX1 = MIN( 2.0_rkind, FPM / SPSIG(IS) )
           AUXH = EXP( -1.0_rkind*(AUX1**4.0_rkind) )
           DO ID = 1, MDC
             IF (SPSIG(IS) .GE. (0.7_rkind*FPM)) THEN
               AUX2 = ( UFRIC(IP) * MAX( 0._rkind , MyCOS(SPDIR(ID)-WINDTH) ) )**4
               SWINA = MAX(0._rkind,AUX * AUX2 * AUXH)
               SSINL(IS,ID) = SWINA / SPSIG(IS)
             END IF
           END DO
         END DO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SIN_EXP( IP, WINDTH, ACLOC, SSINE, DSSINE )
         USE DATAPOOL
         IMPLICIT NONE
!
!     *** the exponential growth term by Komen et al. (1984) ***
!
         INTEGER, INTENT(IN)          :: IP
         REAL(rkind)   , INTENT(IN)   :: WINDTH
         REAL(rkind)   , INTENT(IN)   :: ACLOC(MSC,MDC)
         REAL(rkind)   , INTENT(OUT)  :: SSINE(MSC,MDC), DSSINE(MSC,MDC)

         INTEGER                      :: IS, ID
         REAL(rkind)                  :: AUX1, AUX2, AUX3
         REAL(rkind)                  :: SWINB, CINV, COSDIF

         AUX1 = 0.25_rkind * RHOAW
         AUX2 = 28._rkind * UFRIC(IP)

         DO IS = 1, MSC
           CINV = WK(IS,IP)/SPSIG(IS)
           AUX3 = AUX2 * CINV
           DO ID = 1, MDC
             COSDIF        = MyCOS(SPDIR(ID)-WINDTH)
             SWINB         = AUX1 * ( AUX3  * COSDIF - ONE )
             DSSINE(IS,ID) = MAX( ZERO, SWINB * SPSIG(IS) )
             SSINE(IS,ID)  = DSSINE(IS,ID) * ACLOC(IS,ID)
             DSSINE(IS,ID) = 0.d0
           END DO
         END DO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
