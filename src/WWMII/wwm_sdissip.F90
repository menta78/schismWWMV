      SUBROUTINE SDISSIP (F, FL, IJS, IJL, IG, SL, F1MEAN, XKMEAN)

! ----------------------------------------------------------------------

!**** *SDISSIP* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.

!     S.D.HASSELMANN.
!     MODIFIED TO SHALLOW WATER : G. KOMEN , P. JANSSEN
!     OPTIMIZATION : L. ZAMBRESKY
!     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
!     J. BIDLOT   ECMWF  NOVEMBER 2004  REFORMULATION BASED ON XKMEAN
!                                       AND F1MEAN.
!     P. JANSSEN  ECMWF  JANUARY 2006   ADD BOTTOM-INDUCED DISSIPATION.

!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDISSIP (F, FL, IJS, IJL, SL, F1MEAN, XKMEAN)*
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *IG*  - BLOCK NUMBER
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY
!          *F1MEAN* - MEAN FREQUENCY BASED ON 1st MOMENT.
!          *XKMEAN* - MEAN WAVE NUMBER BASED ON 1st MOMENT.

!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       G.KOMEN, S. HASSELMANN AND K. HASSELMANN, ON THE EXISTENCE
!          OF A FULLY DEVELOPED WINDSEA SPECTRUM, JGR, 1984.

! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FR
!      USE YOWMEAN  , ONLY : EMEAN
!      USE YOWPARAM , ONLY : NANG     ,NFRE
!      USE YOWPCONS , ONLY : G        ,ZPI
!      USE YOWSHAL  , ONLY : DEPTH    ,TFAK     ,INDEP
!      USE YOWSTAT  , ONLY : ISHALLO  ,LBIWBK
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, EMEAN, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, &
     &                      IUSTAR, IALPHA, USTARM, TAUT, STAT, &
     &                      DELUST, DELALP, LBIWBK, DEP,&
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NFRE => MSC, &
     &                      NANG => MDC, &
     &                      INDEP => DEP

! ----------------------------------------------------------------------

      DIMENSION TEMP1(IJS:IJL), SDS(IJS:IJL)

! ----------------------------------------------------------------------

      REAL(rkind),DIMENSION(IJS:IJL,NANG,NFRE) :: F,FL,SL
      REAL(rkind),DIMENSION(IJS:IJL) :: F1MEAN, XKMEAN

      PARAMETER (CDIS = 2.1)
      PARAMETER (DELTA = 0.6)
      PARAMETER (ALPH_B_J = 1.0)
      PARAMETER (GAM_B_J = 0.8)
      PARAMETER (COEF_B_J=2*ALPH_B_J)
      !REAL ZHOOK_HANDLE

      !IF (LHOOK) CALL DR_HOOK('SDISSIP',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
!*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
!        --------------------------------------------------------------

      IF (ISHALLO.EQ.1) THEN
        CONSD = -CDIS*ZPI**9/G**4
        DO IJ=IJS,IJL
          SDS(IJ) = CONSD*F1MEAN(IJ)*EMEAN(IJ)**2*F1MEAN(IJ)**8
        ENDDO
        DO M=1,NFRE
          DO IJ=IJS,IJL
            X         = (FR(M)/F1MEAN(IJ))**2
            TEMP1(IJ) = SDS(IJ)*( (1.-DELTA)*X + DELTA*X**2)
          ENDDO
          DO K=1,NANG
            DO IJ=IJS,IJL
              SL(IJ,K,M) = SL(IJ,K,M)+TEMP1(IJ)*F(IJ,K,M)
              FL(IJ,K,M) = FL(IJ,K,M)+TEMP1(IJ)
            ENDDO
          ENDDO
        ENDDO
      ELSE
!SHALLOW
        CONSS = -CDIS*ZPI
        DO IJ=IJS,IJL
          SDS(IJ) = CONSS*F1MEAN(IJ)*EMEAN(IJ)**2*XKMEAN(IJ)**4
        ENDDO

        DO M=1,NFRE
          DO IJ=IJS,IJL
            !X         = TFAK(INDEP(IJ),M)/XKMEAN(IJ)
            X         = WK(IJ,M)/XKMEAN(IJ)
            TEMP1(IJ) = SDS(IJ)*( (1.-DELTA)*X + DELTA*X**2)
          ENDDO
          DO K=1,NANG
            DO IJ=IJS,IJL
              SL(IJ,K,M) = SL(IJ,K,M)+TEMP1(IJ)*F(IJ,K,M)
              FL(IJ,K,M) = FL(IJ,K,M)+TEMP1(IJ)
            ENDDO
          ENDDO
        ENDDO
!
!*    2. COMPUTATION OF BOTTOM-INDUCED DISSIPATION COEFFICIENT.
!        ----------- -- -------------- -----------------------
!
!        (FOLLOWING BATTJES-JANSSEN AND BEJI)
        DEPTHTRS=50.
        IF(LBIWBK) THEN
          DO IJ=IJS,IJL
             IF(DEPTH(IJ,IG).LT.DEPTHTRS) THEN
               EMAX = (GAM_B_J*DEPTH(IJ,IG))**2/16.
               ALPH = 2.*EMAX/(EMEAN(IJ))
               ARG  = MIN(ALPH,50.)
!!!!!!!! test an iterative scheme
!!!!!!!! if it works we might want to introduce a table
               Q_OLD = EXP(-ARG)
               DO IC=1,15
                 Q = EXP(-ARG*(1.-Q_OLD))
                 REL_ERR=ABS(Q-Q_OLD)/Q_OLD
                 IF(REL_ERR.LT.0.01) EXIT
                 Q_OLD = Q
               ENDDO
               SDS(IJ) = COEF_B_J*ALPH*Q*F1MEAN(IJ)
             ENDIF
          ENDDO 
      
          DO M=1,NFRE
             DO K=1,NANG
                DO IJ=IJS,IJL
                  IF(DEPTH(IJ,IG).LT.DEPTHTRS) THEN
                    SL(IJ,K,M) = SL(IJ,K,M)-SDS(IJ)*F(IJ,K,M)
                    FL(IJ,K,M) = FL(IJ,K,M)-SDS(IJ)
                  ENDIF
                ENDDO
             ENDDO
          ENDDO
        ENDIF
     
!SHALLOW
      ENDIF

      IF (LHOOK) CALL DR_HOOK('SDISSIP',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE SDISSIP
! ----------------------------------------------------------------------

      SUBROUTINE SDISSIP_WWM (IP, F, FL, IG, SL, EMEAN, F1MEAN, XKMEAN)

! ----------------------------------------------------------------------

!**** *SDISSIP* - COMPUTATION OF DISSIPATION SOURCE FUNCTION.

!     S.D.HASSELMANN.
!     MODIFIED TO SHALLOW WATER : G. KOMEN , P. JANSSEN
!     OPTIMIZATION : L. ZAMBRESKY
!     J. BIDLOT   ECMWF  FEBRUARY 1997   ADD SL IN SUBROUTINE CALL
!     J. BIDLOT   ECMWF  NOVEMBER 2004  REFORMULATION BASED ON XKMEAN
!                                       AND F1MEAN.
!     P. JANSSEN  ECMWF  JANUARY 2006   ADD BOTTOM-INDUCED DISSIPATION.

!*    PURPOSE.
!     --------
!       COMPUTE DISSIPATION SOURCE FUNCTION AND STORE ADDITIVELY INTO
!       NET SOURCE FUNCTION ARRAY. ALSO COMPUTE FUNCTIONAL DERIVATIVE
!       OF DISSIPATION SOURCE FUNCTION.

!**   INTERFACE.
!     ----------

!       *CALL* *SDISSIP (F, FL, IJS, IJL, SL, F1MEAN, XKMEAN)*
!          *F*   - SPECTRUM.
!          *FL*  - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE
!          *IJS* - INDEX OF FIRST GRIDPOINT
!          *IJL* - INDEX OF LAST GRIDPOINT
!          *IG*  - BLOCK NUMBER
!          *SL*  - TOTAL SOURCE FUNCTION ARRAY
!          *F1MEAN* - MEAN FREQUENCY BASED ON 1st MOMENT.
!          *XKMEAN* - MEAN WAVE NUMBER BASED ON 1st MOMENT.

!     METHOD.
!     -------

!       SEE REFERENCES.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       G.KOMEN, S. HASSELMANN AND K. HASSELMANN, ON THE EXISTENCE
!          OF A FULLY DEVELOPED WINDSEA SPECTRUM, JGR, 1984.

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FR
!      USE YOWMEAN  , ONLY : EMEAN
!      USE YOWPARAM , ONLY : NANG     ,NFRE
!      USE YOWPCONS , ONLY : G        ,ZPI
!      USE YOWSHAL  , ONLY : DEPTH    ,TFAK     ,INDEP
!      USE YOWSTAT  , ONLY : ISHALLO  ,LBIWBK
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!      USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPRM

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, &
     &                      IUSTAR, IALPHA, USTARM, TAUT, STAT, &
     &                      DELUST, DELALP, LBIWBK, DEP,&
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NFRE => MSC, &
     &                      NANG => MDC, &
     &                      INDEP => DEP

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IG, IP

      REAL(rkind) :: TEMP1, SDS

! ----------------------------------------------------------------------

      REAL(rkind) :: F(NANG,NFRE)
      REAL(rkind) :: FL(NANG,NFRE)
      REAL(rkind) :: SL(NANG,NFRE)
      REAL(rkind) :: F1MEAN, XKMEAN, EMEAN

      REAL(rkind), PARAMETER :: CDIS = 2.1
      REAL(rkind), PARAMETER :: DELTA = 0.6
      REAL(rkind), PARAMETER :: ALPH_B_J = 1.0
      REAL(rkind), PARAMETER :: GAM_B_J = 0.8
      REAL(rkind), PARAMETER :: COEF_B_J=2*ALPH_B_J

      REAL(rkind) :: CONSD, CONSS, X, DEPTHTRS, Q, ARG, EMAX, ALPH, Q_OLD, REL_ERR
      INTEGER     :: IJ, I, J, K, L, M, N, IC
      !REAL(KIND=JPRB) ::  ZHOOK_HANDLE

      !IF (LHOOK) CALL DR_HOOK('SDISSIP',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    1. ADDING DISSIPATION AND ITS FUNCTIONAL DERIVATIVE TO NET SOURCE
!*       FUNCTION AND NET SOURCE FUNCTION DERIVATIVE.
!        --------------------------------------------------------------

      IF (ISHALLO.EQ.1) THEN
        CONSD = -CDIS*ZPI**9/G**4
        SDS = CONSD*F1MEAN*EMEAN**2*F1MEAN**8
        DO M=1,NFRE
          X         = (FR(M)/F1MEAN)**2
          TEMP1 = SDS*( (1.-DELTA)*X + DELTA*X**2)
          DO K=1,NANG
            SL(K,M) = SL(K,M)+TEMP1*F(K,M)
            FL(K,M) = FL(K,M)+TEMP1
          ENDDO
        ENDDO
      ELSE
!SHALLOW
        CONSS = -CDIS*ZPI
        SDS   = CONSS*F1MEAN*EMEAN**2*XKMEAN**4

        DO M=1,NFRE
!            X         = TFAK(INDEP(IJ),M)/XKMEAN(IJ)
          X         = WK(IP,M)/XKMEAN
          TEMP1     = SDS*( (1.-DELTA)*X + DELTA*X**2)
!          WRITE(STAT%FHNDL,'(I10,5D15.5)') M, X, XKMEAN, TEMP1, WK(IP,M)
          DO K=1,NANG
            SL(K,M) = SL(K,M)+TEMP1*F(K,M)
            FL(K,M) = FL(K,M)+TEMP1
            !write(DBG%FHNDL,'(2I10,10F15.8)') M,K,FL(K,M),TEMP1,4*SQRT(EMEAN)
          ENDDO
        ENDDO
!
!*    2. COMPUTATION OF BOTTOM-INDUCED DISSIPATION COEFFICIENT.
!        ----------- -- -------------- -----------------------
!
        DEPTHTRS=50.
        IF(LBIWBK .or. .false.) THEN
           IF(DEP(IP).LT.DEPTHTRS) THEN
             EMAX = (GAM_B_J*DEP(IP))**2/16.
             ALPH = 2.*EMAX/(EMEAN)
             ARG  = MIN(ALPH,50._rkind)
!!!!!!!! test an iterative scheme
!!!!!!!! if it works we might want to introduce a table
             Q_OLD = EXP(-ARG)
             DO IC=1,15
               Q = EXP(-ARG*(1.-Q_OLD))
               REL_ERR=ABS(Q-Q_OLD)/Q_OLD
               IF(REL_ERR.LT.0.01) EXIT
               Q_OLD = Q
             ENDDO
             SDS = COEF_B_J*ALPH*Q*F1MEAN
           ENDIF

          DO M=1,NFRE
             DO K=1,NANG
                IF(DEP(IP).LT.DEPTHTRS) THEN
                  SL(K,M) = SL(K,M)-SDS*F(K,M)
                  FL(K,M) = FL(K,M)-SDS
                ENDIF
             ENDDO
          ENDDO
        ENDIF

      ENDIF

!      IF (LHOOK) CALL DR_HOOK('SDISSIP',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE SDISSIP_WWM
