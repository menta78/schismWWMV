#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ST4_PRE (IP, ACLOC, IMATRA, IMATDA, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
      USE DATAPOOL
      USE W3SRC4MD
      IMPLICIT NONE

      INTEGER, INTENT(IN)        :: IP
      REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)

      REAL(rkind), INTENT(OUT)   :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
      REAL(rkind), INTENT(OUT)   :: SSINE(MSC,MDC),DSSINE(MSC,MDC) 
      REAL(rkind), INTENT(OUT)   :: SSDS(MSC,MDC),DSSDS(MSC,MDC)
      REAL(rkind), INTENT(OUT)   :: SSNL4(MSC,MDC),DSSNL4(MSC,MDC)
      REAL(rkind), INTENT(OUT)   :: SSINL(MSC,MDC)

      INTEGER      :: IS, ID, ITH, IK, IS0

      REAL(rkind)  :: AWW3(NSPEC), IMATRAWW3(MSC,MDC), IMATDAWW3(MSC,MDC)
      REAL(rkind)  :: IMATDA1D(NSPEC), IMATRA1D(NSPEC), BRLAMBDA(NSPEC)
      REAL(rkind)  :: WN2(MSC*MDC), WHITECAP(1:4)

      REAL(rkind)  :: ETOT, FAVG, FMEAN1, WNMEAN, AS, SUMACLOC, FAVGWS
      REAL(rkind)  :: TAUWAX, TAUWAY, AMAX, FPM, WIND10, WINDTH
      REAL(rkind)  :: HS,SME01,SME10,KME01,KMWAM,KMWAM2

      IMATRA = 0.d0 
      IMATDA = 0.d0

      DO IS = 1, MSC
        DO ID = 1, MDC
          AWW3(ID + (IS-1) * MDC) = ACLOC(IS,ID) * CG(IS,IP)
        END DO
      END DO

      DO IK=1, NK
        WN2(1+(IK-1)*NTH) = WK(IK,IP)
      END DO

      DO IK=1, NK
        IS0    = (IK-1)*NTH
        DO ITH=2, NTH
          WN2(ITH+IS0) = WN2(1+IS0)
        END DO
      END DO
!
! wind input
!
      TAUWX(IP)  = ZERO
      TAUWY(IP)  = ZERO               
      SSINL      = ZERO
      MSC_HF(IP) = MSC
      AS         = 0. 

      IF (MESIN .GT. 0) THEN

        CALL SET_WIND( IP, WIND10, WINDTH )
        CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
        CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)
        CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), ETOT, FMEAN(IP), FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS(IP))
        IF (ETOT .LT. THR .AND. WIND10 .GT. THR) CALL SIN_LIN_CAV(IP,WINDTH,FPM,SSINL)
        CALL W3SIN4 ( IP, AWW3, CG(:,IP), WN2, WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, IMATRA1D, IMATDA1D, LLWS, BRLAMBDA)
        CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), ETOT, FMEAN(IP), FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS(IP))  
        CALL W3SIN4 ( IP, AWW3, CG(:,IP), WN2, WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, IMATRA1D, IMATDA1D, LLWS, BRLAMBDA)
#ifdef DEBUG
        WRITE(740+myrank,*) 'min/max/sum(VSIN)=', minval(IMATRA1D), maxval(IMATRA1D), sum(IMATRA1D)
        WRITE(740+myrank,*) 'min/max/sum(VDIN)=', minval(IMATDA1D), maxval(IMATDA1D), sum(IMATDA1D)
#endif
        CALL ONED2TWOD(IMATRA1D,IMATRAWW3)
        CALL ONED2TWOD(IMATDA1D,IMATDAWW3)
        DO ID = 1, MDC
          SSINE(:,ID)  = IMATRAWW3(:,ID) / CG(:,IP) 
          DSSINE(:,ID) = IMATDAWW3(:,ID)
        END DO

      ENDIF

      IF (MESNL .GT. 0) CALL SNL41(IP, KMWAM, ACLOC, IMATRA, IMATDA, SSNL4, DSSNL4)

      IF (MESDS .GT. 0) THEN
        CALL W3SDS4(AWW3,WK(:,IP),CG(:,IP),UFRIC(IP),USTDIR(IP),DEP(IP),IMATRA1D,IMATDA1D,BRLAMBDA,WHITECAP)
        CALL ONED2TWOD(IMATRA1D,IMATRAWW3)
        CALL ONED2TWOD(IMATDA1D,IMATDAWW3)
        DO ID = 1, MDC
          SSDS(:,ID)  = IMATRAWW3(:,ID) / CG(:,IP)
          DSSDS(:,ID) = IMATDAWW3(:,ID) 
        END DO
      ENDIF

      IMATRA = SSINL + SSINE + SSNL4 + SSDS
      IMATDA = DSSINE + DSSNL4 + DSSDS
!
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ST4_POST (IP, ACLOC, SSINE, DSSINE, SSDS, DSSDS, SSINL)
        USE DATAPOOL
        USE W3SRC4MD
        IMPLICIT NONE
       
        INTEGER, INTENT(IN)        :: IP
        REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)
       
        REAL(rkind), INTENT(OUT)   :: SSINE(MSC,MDC),DSSINE(MSC,MDC)
        REAL(rkind), INTENT(OUT)   :: SSDS(MSC,MDC),DSSDS(MSC,MDC)
        REAL(rkind), INTENT(OUT)   :: SSINL(MSC,MDC)

        INTEGER                    :: IS, ID, IK, ITH, ITH2, IS0

        REAL(rkind)                :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
        REAL(rkind)                :: AWW3(NSPEC), WN2(MSC*MDC), BRLAMBDA(NSPEC)
        REAL(rkind)                :: IMATDA1D(NSPEC), IMATRA1D(NSPEC), TMP_DS(MSC)
        REAL(rkind)                :: IMATRAWW3(MSC,MDC), IMATDAWW3(MSC,MDC)

        REAL(rkind)                :: ETOT, FAVG, FMEAN1, WNMEAN, AS, FAVGWS
        REAL(rkind)                :: TAUWAX, TAUWAY, AMAX, WIND10, WINDTH
        REAL(rkind)                :: WHITECAP(1:4), SUMACLOC, FPM

        DO IS = 1, MSC
          DO ID = 1, MDC
            AWW3(ID + (IS-1) * MDC) = ACLOC(IS,ID) * CG(IS,IP)
          END DO
        END DO

        DO IK=1, NK
          WN2(1+(IK-1)*NTH) = WK(IK,IP)
        END DO
        DO IK=1, NK
          IS0    = (IK-1)*NTH
          DO ITH=2, NTH
            WN2(ITH+IS0) = WN2(1+IS0)
          END DO
        END DO

! wind input
        AS      = 0.
        MSC_HF(IP) = MSC
        CALL SET_WIND( IP, WIND10, WINDTH )
        CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
        CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), ETOT, FMEAN(IP), FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS(IP))
        IF (ETOT .LT. THR .AND. WIND10 .GT. THR) CALL SIN_LIN_CAV(IP,WINDTH,FPM,IMATRA,SSINL)
        CALL W3SIN4 ( IP, AWW3, CG(:,IP), WN2,  WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, IMATRA1D, IMATDA1D, LLWS, BRLAMBDA)
        CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), ETOT, FMEAN(IP), FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS(IP))
        CALL ONED2TWOD(IMATRA1D,IMATRAWW3)
        CALL ONED2TWOD(IMATDA1D,IMATDAWW3)
        DO ID = 1, MDC
          SSINE(:,ID)  = IMATRAWW3(:,ID) / CG(:,IP)
          DSSINE(:,ID) = IMATDAWW3(:,ID)
        END DO

! dissipation 
        CALL W3SDS4(AWW3,WK(:,IP),CG(:,IP),UFRIC(IP),USTDIR(IP),DEP(IP),IMATRA1D,IMATDA1D,BRLAMBDA,WHITECAP)
        CALL ONED2TWOD(IMATRA1D,IMATRAWW3)
        CALL ONED2TWOD(IMATDA1D,IMATDAWW3)
        DO ID = 1, MDC
          SSDS(:,ID)   = IMATRAWW3(:,ID) / CG(:,IP)
          DSSDS(:,ID)  = IMATDAWW3(:,ID)
        END DO
        IMATRA = IMATRA + SSDS
        IMATDA = IMATDA + DSSDS

! missing high freq. tail contribution -> 2do

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
