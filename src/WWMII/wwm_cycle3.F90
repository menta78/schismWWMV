#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CYCLE3 (IP, ACLOC, IMATRA, IMATDA)
         USE DATAPOOL

         INTEGER, INTENT(IN)        :: IP

         REAL(rkind), INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)

         CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2) 

         IF (MESIN .GT. 0) THEN
           CALL SET_WIND( IP, WIND10, WINDTH )
           CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
           IF (.NOT. LINID) CALL SIN_LIN_CAV( IP, WINDTH, FPM, SSINL)
           CALL SIN_EXP( IP, WINDTH, ACLOC, SSINE )
         ENDIF

         IF (MESDS .GT. 0) CALL SDSCYCLE3 ( IP, KMWAM, SME10, ETOT, ACLOC, SSDS )
         IF (MESNL .GT. 0) CALL SNL4_NEW  (IP, KMWAM, ACLOC, SSNL4, DSSNL4)
         IF (MESTR .GT. 0 .AND. ISHALLOW(IP) .EQ. 1) CALL TRIADSWAN_NEW (IP,HS,SME01,ACLOC,SSNL3)
         IF (MESBF .GT. 0 .AND. ISHALLOW(IP) .EQ. 1) CALL SDS_SWB(IP,SME01,KMWAM,ETOT,HS,ACLOC,SSBR)
         IF (MESBF .GT. 0 .AND. ISHALLOW(IP) .EQ. 1) CALL SDS_BOTF(IP,ACLOC,SSBF)

         IF (ICOMP .LT. 2) THEN
           IMATRA = SSINL + SSDS + SSNL4 + SSNL3 + SSBR + SSBF
           IMATDA =        DSSDS + DSSNL4 + DSSNL3 + DSSBR + DSSBF
         ELSE
           IMATRA = MIN(ZERO, SSINL + SSDS + SSNL4 + SSNL3 + SSBR + SSBF) ! Patankar Rules ...
           IMATDA = MAX(ZERO, DSSDS + DSSNL4 + DSSNL3 + DSSBR + DSSBF)
         ENDIF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SDSCYCLE3( IP, KMESPC, SMESPC, ETOT, ACLOC, SSDS, DSSDS )
!
!     Cycle 3 dissipation 
!
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)           :: IP
         REAL(rkind)   , INTENT(IN)    :: KMESPC, SMESPC, ETOT
         REAL(rkind)   , INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind)   , INTENT(OUT)   :: SSDS(MSC,MDC), DSSDS(MSC,MDC)

         INTEGER       :: IS, ID

         REAL(rkind)    :: CDS, ALPHA_PM, FAC
         REAL(rkind)    :: STP_OV, STP_PM, N2
!
         ALPHA_PM  =  3.02E-3
         CDS       = -2.36E-5

         STP_OV = KMESPC * SQRT(ETOT)
         STP_PM = SQRT(ALPHA_PM)

         N2     = 4

         FAC    = CDS * (STP_OV / STP_PM)**N2

         DO IS = 1, MSC
           DSSDS(IS,:) = FAC * SMESPC * (WK(IP,IS)/KMESPC)
           SSDS(IS,:)  = DSSDS(IS,:) * ACLOC(IS,:)
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
         REAL(rkind)                  :: SWINB, CINV, COSDIF, SFIE(MSC,MDC)

         AUX1 = 0.25_rkind * RHOAW
         AUX2 = 28._rkind * UFRIC(IP)

         DO IS = 1, MSC
           CINV = WK(IP,IS)/SPSIG(IS)
           AUX3 = AUX2 * CINV
           DO ID = 1, MDC
             COSDIF = MyCOS(SPDIR(ID)-WINDTH)
             SWINB = AUX1 * ( AUX3  * COSDIF - ONE )
             DSSINE(IS,ID) = MAX( ZERO, SWINB * SPSIG(IS) )
             SSINE(IS,ID) = DSSINE(IS,ID) * ACLOC(IS,ID)
           END DO
         END DO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
