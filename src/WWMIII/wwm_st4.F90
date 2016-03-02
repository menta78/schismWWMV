#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ST4_PRE (IP, ACLOC, IMATRA, IMATDA)
      USE DATAPOOL
      USE W3SRC4MD
      IMPLICIT NONE

      INTEGER, INTENT(IN)        :: IP

      REAL(rkind), INTENT(OUT)   :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
      REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)

      INTEGER      :: IS, ID, ITH, IK, IS0
     
      REAL(rkind)  :: SSINE(MSC,MDC),DSSINE(MSC,MDC) 
      REAL(rkind)  :: SSDS(MSC,MDC),DSSDS(MSC,MDC)
      REAL(rkind)  :: SSNL4(MSC,MDC),DSSNL4(MSC,MDC)
      REAL(rkind)  :: SSINL(MSC,MDC)

      REAL(rkind)  :: AWW3(NSPEC), IMATRAWW3(MSC,MDC), IMATDAWW3(MSC,MDC)
      REAL(rkind)  :: IMATDA1D(NSPEC), IMATRA1D(NSPEC), BRLAMBDA(NSPEC)
      REAL(rkind)  :: WN2(MSC*MDC), WHITECAP(1:4)

      REAL(rkind)  :: XRR, XPP, XFLT, XREL, FACP, XFILT, FAVGWS
      REAL(rkind)  :: ETOT, FAVG, FMEAN1, WNMEAN, AS, SUMACLOC
      REAL(rkind)  :: TAUWAX, TAUWAY, AMAX, FPM, WIND10, WINDTH, KMWAM

      XPP     = 0.15
      XRR     = 0.10
      XFILT  = 0.05
      XPP     = MAX ( 1.E-6_rkind , XPP )
      XRR     = MAX ( 1.E-6_rkind , XRR )
      XREL   = XRR
      XFILT  = MAX ( ZERO , XFILT )
      XFLT   = XFILT
      FACP   = 2*XPP / PI2 * 0.62E-3 * PI2**4 / G9**2

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

      CALL SET_WIND( IP, WIND10, WINDTH )

      TAUWX(IP) = ZERO
      TAUWY(IP) = ZERO               

      SUMACLOC = SUM(ACLOC)

      IF (IOBP(IP) .EQ. 0) THEN
        IF (SUMACLOC .LT. THR .AND. WIND10 .GT. THR .AND. .NOT. LINID) THEN
          CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
          CALL SIN_LIN_CAV(IP,WINDTH,FPM,IMATRA,SSINL)
        ELSE
          MSC_HF(IP) = MSC
          AS      = 0. 
          CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), ETOT, FAVG, FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FAVGWS)
          CALL W3SIN4 ( IP, AWW3, CG(:,IP), WN2,  WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, IMATRA1D, IMATDA1D, LLWS, BRLAMBDA)
          CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), ETOT, FAVG, FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FAVGWS)  
          CALL W3SIN4 ( IP, AWW3, CG(:,IP), WN2, WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, IMATRA1D, IMATDA1D, LLWS, BRLAMBDA) 
          CALL ONED2TWOD(IMATRA1D,IMATRAWW3)
          CALL ONED2TWOD(IMATDA1D,IMATDAWW3)
          DO ID = 1, MDC
            SSINE(:,ID) = IMATRAWW3(:,ID) / CG(:,IP) 
            DSSINE(:,ID) = IMATDAWW3(:,ID)
          END DO
          IMATRA = IMATRA + SSINE
          IMATDA = IMATDA + DSSINE  
        END IF
        CALL SNL41(IP, KMWAM, ACLOC, IMATRA, IMATDA, SSNL4, DSSNL4)
        IMATRA = IMATRA + SSNL4 
        IMATDA = IMATDA + DSSNL4 
        CALL W3SDS4(AWW3,WK(:,IP),CG(:,IP),UFRIC(IP),USTDIR(IP),DEP(IP),IMATRA1D,IMATDA1D,BRLAMBDA,WHITECAP)
        CALL ONED2TWOD(IMATRA1D,IMATRAWW3)
        CALL ONED2TWOD(IMATDA1D,IMATDAWW3)
        DO ID = 1, MDC
          SSDS(:,ID)   = IMATRAWW3(:,ID) / CG(:,IP)
          DSSDS(:,ID)  = IMATDAWW3(:,ID) 
        END DO
        IMATRA = IMATRA + SSDS 
        IMATDA = IMATDA + DSSDS 
      ELSE
        IF (LSOUBOUND) THEN
          IF (IOBP(IP) .NE. 2) THEN
            IF (SUMACLOC .LT. THR .AND. WIND10 .GT. THR .AND. .NOT. LINID) THEN
              CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
              CALL SIN_LIN_CAV(IP,WINDTH,FPM,IMATRA,SSINL)
            ELSE
              MSC_HF(IP) = MSC
              AS      = 0.
              CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), ETOT, FAVG, FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FAVGWS)
              CALL W3SIN4 ( IP, AWW3, CG(:,IP), WN2,  WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, IMATRA1D, IMATDA1D, LLWS, BRLAMBDA)
              CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), ETOT, FAVG, FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FAVGWS)
              CALL W3SIN4 ( IP, AWW3, CG(:,IP), WN2, WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, IMATRA1D, IMATDA1D, LLWS, BRLAMBDA)
              CALL ONED2TWOD(IMATRA1D,IMATRAWW3)
              CALL ONED2TWOD(IMATDA1D,IMATDAWW3)
              DO ID = 1, MDC
                SSINE(:,ID) = IMATRAWW3(:,ID) / CG(:,IP)
                DSSINE(:,ID) = IMATDAWW3(:,ID)
              END DO
              IMATRA = IMATRA + SSINE
              IMATDA = IMATDA + DSSINE
            END IF
            CALL SNL41(IP, KMWAM, ACLOC, IMATRA, IMATDA, SSNL4, DSSNL4)
            IMATRA = IMATRA + SSNL4
            IMATDA = IMATDA + DSSNL4
            CALL W3SDS4(AWW3,WK(:,IP),CG(:,IP),UFRIC(IP),USTDIR(IP),DEP(IP),IMATRA1D,IMATDA1D,BRLAMBDA,WHITECAP)
            CALL ONED2TWOD(IMATRA1D,IMATRAWW3)
            CALL ONED2TWOD(IMATDA1D,IMATDAWW3)
            DO ID = 1, MDC
              SSDS(:,ID)   = IMATRAWW3(:,ID) / CG(:,IP)
              DSSDS(:,ID)  = IMATDAWW3(:,ID)
            END DO
            IMATRA = IMATRA + SSDS
            IMATDA = IMATDA + DSSDS
          ENDIF
        ENDIF
      ENDIF

      IMATDA = 0.d0

      IF (IP == TESTNODE) THEN
        WRITE(*,'(A20,6E20.10)') 'LINEAR INPUT', SUM(SSINL), MINVAL(SSINL), MAXVAL(SSINL)
        WRITE(*,'(A20,6E20.10)') 'WAVE ACTION', SUM(ACLOC), MINVAL(ACLOC), MAXVAL(ACLOC)
        WRITE(*,'(A20,6E20.10)') 'EXP INPUT', SUM(SSINE), SUM(DSSINE), MINVAL(SSINE), MAXVAL(SSINE), MINVAL(DSSINE), MAXVAL(DSSINE)
        WRITE(*,'(A20,6E20.10)') 'WHITECAP', SUM(SSDS), SUM(DSSDS), MINVAL(SSDS), MAXVAL(SSDS), MINVAL(DSSDS), MAXVAL(DSSDS)
        WRITE(*,'(A20,6E20.10)') 'SNL4', SUM(SSNL4), SUM(DSSNL4), MINVAL(SSNL4), MAXVAL(SSNL4), MINVAL(DSSNL4), MAXVAL(DSSNL4)
        WRITE(*,'(A20,6E20.10)') 'TOTAL SOURCE TERMS', SUM(IMATRA), SUM(IMATDA), MINVAL(IMATRA), MAXVAL(IMATRA), MINVAL(IMATDA), MAXVAL(IMATDA)
      ENDIF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
