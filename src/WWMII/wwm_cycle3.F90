#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CYCLE3 (IP, ACLOC, IMATRA, IMATDA)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IP

         REAL(rkind), INTENT(OUT)   :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)

         INTEGER                    :: IS, ID

         REAL(rkind)                :: NEWAC(MSC,MDC), SSINL(MSC,MDC)
         REAL(rkind)                :: SSINE(MSC,MDC),DSSINE(MSC,MDC)
         REAL(rkind)                :: SSDS(MSC,MDC),DSSDS(MSC,MDC)
         REAL(rkind)                :: SSNL4(MSC,MDC),DSSNL4(MSC,MDC)
         REAL(rkind)                :: SSNL3(MSC,MDC),DSSNL3(MSC,MDC)
         REAL(rkind)                :: SSBR(MSC,MDC),DSSBR(MSC,MDC)
         REAL(rkind)                :: SSBF(MSC,MDC),DSSBF(MSC,MDC)
         REAL(rkind)                :: SSBRL(MSC,MDC),DSSBRL(MSC,MDC)
         REAL(rkind)                :: ETOT,SME01,SME10,KME01,KMWAM,KMWAM2,HS,WIND10
         REAL(rkind)                :: ETAIL,EFTAIL,EMAX,LIMAC,NEWDAC,MAXDAC,FPM,WINDTH
         REAL(rkind)                :: RATIO,LIMFAC

         CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2) 

         IF (MESIN .GT. 0) THEN
           CALL SET_WIND( IP, WIND10, WINDTH )
           CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
           IF (.NOT. LINID) CALL SIN_LIN( IP, WINDTH, FPM, SSINL)
           CALL SIN_EXP( IP, WINDTH, ACLOC, SSINE, DSSINE )
         ENDIF

         IF (MESDS .GT. 0) CALL SDS_CYCLE3_NEW ( IP, KMWAM, SME10, ETOT, ACLOC, SSDS, DSSDS )
         IF (MESNL .GT. 0) CALL SNL4_NEW  (IP, KMWAM, ACLOC, SSNL4, DSSNL4)

         IF (ISHALLOW(IP) .EQ. 1) THEN
           !IF (MESTR .GT. 0) CALL TRIADSWAN_NEW2 (IP,HS,SME01,ACLOC,SSNL3, DSSNL3)
           !IF (MESBF .GT. 0) CALL SDS_SWB_NEW(IP,SME01,KMWAM,ETOT,HS,ACLOC,SSBR,DSSBR)
           !IF (MESBF .GT. 0) CALL SDS_BOTF_NEW(IP,ACLOC,SSBF,DSSBF)
         ENDIF

         IMATRA = SSINL + SSINE +  SSDS +  SSNL4 +  SSNL3 
         IMATDA =        DSSINE + DSSDS + DSSNL4 + DSSNL3 

         RETURN

         DO IS = 1, MSC
           MAXDAC   = LIMFAK*0.0081_rkind/(TWO*SPSIG(IS)*WK(IP,IS)**3*CG(IP,IS))
           DO ID = 1, MDC
             NEWDAC        = IMATRA(IS,ID)*DT4A/MAX((ONE-DT4A*IMATDA(IS,ID)),ONE) 
             IMATRA(IS,ID) = MIN(ABS(NEWDAC),MAXDAC)/DT4A ! This is now the source term ... right hand side
             LIMFAC        = MIN(ONE,ABS(SIGN(MAXDAC/DT4A,NEWDAC/DT4A))/MAX(THR,ABS(IMATRA(IS,ID))))
             IMATDA(IS,ID) = LIMFAC * IMATDA(IS,ID) ! This is the new source term ... diagonal part 
           ENDDO
         ENDDO

         IMATRA = MAX(ZERO, IMATRA +  SSBR +  SSBF) 
         IMATDA = MIN(ZERO, IMATDA + DSSBR + DSSBF)

         NEWAC = ACLOC + IMATRA*DT4A/MAX((ONE-DT4A*IMATDA),ONE)
         ETOT   = ZERO 
         EFTAIL = ONE / (PTAIL(1)-ONE)
         HS = 4._rkind*SQRT(ETOT)
         EMAX = 1._rkind/16._rkind * (HMAX(IP))**2 ! HMAX is defined in the breaking routine or has some default value
         IF (ETOT .GT. EMAX) THEN
           RATIO  = EMAX/ETOT
           SSBRL  = ACLOC*(RATIO-ONE)/DT4A
           DSSBRL = (RATIO-ONE)/DT4A 
         END IF

         IMATRA = MAX(ZERO, IMATRA +  SSBRL) 
         IMATDA = MIN(ZERO, IMATDA + DSSBRL)

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
