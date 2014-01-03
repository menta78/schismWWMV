      SUBROUTINE FEMEANWS (F, IJS, IJL, USNEW, THWNEW, EM, FM, XLLWS)

! ----------------------------------------------------------------------

!**** *FEMEANWS* - COMPUTATION OF MEAN ENERGY, MEAN FREQUENCY 
!                  FOR WINDSEA PART OF THE SPECTRUM AS DETERMINED
!                  BY THE EMPIRICAL LAW BASED ON WAVE AGE AND
!                  THE DIRECTIOn WITH RESPECT TO THE WIND DIRECTION
!                  (SEE LLWS)

!*    PURPOSE.
!     --------

!       COMPUTE MEAN FREQUENCY AT EACH GRID POINT FOR PART OF THE
!       SPECTRUM WHERE LLWS IS TRUE OR THE WINDSEA PARAMETRIC LAW
!       APPLIES.

!**   INTERFACE.
!     ----------

!       *CALL* *FEMEANWS (F, IJS, IJL, EM, FM)*
!              *F*      - SPECTRUM.
!              *IJS*    - INDEX OF FIRST GRIDPOINT
!              *IJL*    - INDEX OF LAST GRIDPOINT
!              *USNEW*  - FRICTION VELOCITY
!              *THWNEW* - WIND DIRECTION
!              *EM*     - MEAN WAVE ENERGY (OUTPUT)
!              *FM*     - MEAN WAVE FREQUENCY (OUTPUT)

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMOFR  ,DELTH    ,
!     &                WETAIL    ,FRTAIL     ,TH    ,C     ,FRIC    
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!      USE YOWPARAM , ONLY : NANG     ,NFRE
!      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN
       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, &
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP


! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IJ,M,K,IJS,IJL
      REAL(rkind) :: DELT25, DELT2, CM, CHECKTA
      REAL(rkind) :: F(IJS:IJL,NANG,NFRE)
      REAL(rkind),DIMENSION(IJS:IJL) :: THWNEW,USNEW
      REAL(rkind),DIMENSION(IJS:IJL) :: TEMP2, EM, FM, THRESHOLD
      REAL(rkind), DIMENSION(IJS:IJL,NANG,NFRE) :: XLLWS

      !REAL ZHOOK_HANDLE

      !IF (LHOOK) CALL DR_HOOK('FEMEANWS',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      DO IJ=IJS,IJL
        EM(IJ) = EPSMIN
        FM(IJ) = EPSMIN
      ENDDO

      DELT25 = WETAIL*FR(NFRE)*DELTH
      DELT2 = FRTAIL*DELTH


!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------
      
      DO M=1,NFRE
        K = 1
        DO IJ =IJS,IJL
           TEMP2(IJ) = F(IJ,K,M)*XLLWS(IJ,K,M)
        ENDDO
        DO K=2,NANG
          DO IJ=IJS,IJL
            TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)*XLLWS(IJ,K,M)
          ENDDO
        ENDDO
        DO IJ=IJS,IJL
          EM(IJ) = EM(IJ)+TEMP2(IJ)*DFIM(M)
          FM(IJ) = FM(IJ)+DFIMOFR(M)*TEMP2(IJ)
        ENDDO
      ENDDO

!*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
!*       NORMALIZE WITH TOTAL ENERGY.
!        ------------------------------------------

      DO IJ=IJS,IJL
        EM(IJ) = EM(IJ)+DELT25*TEMP2(IJ)
        FM(IJ) = FM(IJ)+DELT2*TEMP2(IJ)
        FM(IJ) = EM(IJ)/FM(IJ)
      ENDDO

      !IF (LHOOK) CALL DR_HOOK('FEMEANWS',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE FEMEANWS
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
      SUBROUTINE FEMEANWS_WWM (IP, F, USNEW, THWNEW, EM, FM, XLLWS)

! ----------------------------------------------------------------------

!**** *FEMEANWS* - COMPUTATION OF MEAN ENERGY, MEAN FREQUENCY
!                  FOR WINDSEA PART OF THE SPECTRUM AS DETERMINED
!                  BY THE EMPIRICAL LAW BASED ON WAVE AGE AND
!                  THE DIRECTIOn WITH RESPECT TO THE WIND DIRECTION
!                  (SEE LLWS)

!*    PURPOSE.
!     --------

!       COMPUTE MEAN FREQUENCY AT EACH GRID POINT FOR PART OF THE
!       SPECTRUM WHERE LLWS IS TRUE OR THE WINDSEA PARAMETRIC LAW
!       APPLIES.

!**   INTERFACE.
!     ----------

!       *CALL* *FEMEANWS (F, IJS, IJL, EM, FM)*
!              *F*      - SPECTRUM.
!              *IJS*    - INDEX OF FIRST GRIDPOINT
!              *IJL*    - INDEX OF LAST GRIDPOINT
!              *USNEW*  - FRICTION VELOCITY
!              *THWNEW* - WIND DIRECTION
!              *EM*     - MEAN WAVE ENERGY (OUTPUT)
!              *FM*     - MEAN WAVE FREQUENCY (OUTPUT)

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMOFR  ,DELTH    ,
!     &                WETAIL    ,FRTAIL     ,TH    ,C     ,FRIC
!      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!      USE YOWPARAM , ONLY : NANG     ,NFRE
!      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN
!      USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPRM

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, &
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IP

      INTEGER :: IJ,M,K
      REAL(rkind) :: DELT25, DELT2, CM, CHECKTA
      REAL(rkind) :: F(NANG,NFRE)
      REAL(rkind) :: THWNEW,USNEW
      REAL(rkind) :: TEMP2, EM, FM, THRESHOLD
      REAL(rkind), DIMENSION(NANG,NFRE) :: XLLWS

      !REAL(KIND=JPRB) :: ZHOOK_HANDLE

      !IF (LHOOK) CALL DR_HOOK('FEMEANWS',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      EM = EPSMIN
      FM = EPSMIN

      DELT25 = WETAIL*FR(NFRE)*DELTH
      DELT2 = FRTAIL*DELTH

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      DO M=1,NFRE
        K = 1
        TEMP2 = F(K,M)*XLLWS(K,M)
        DO K=2,NANG
          TEMP2 = TEMP2+F(K,M)*XLLWS(K,M)
        ENDDO
        EM = EM+TEMP2*DFIM(M)
        FM = FM+DFIMOFR(M)*TEMP2
      ENDDO

!*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
!*       NORMALIZE WITH TOTAL ENERGY.
!        ------------------------------------------
      EM = EM+DELT25*TEMP2
      FM = FM+DELT2*TEMP2
      FM = EM/FM

      !IF (LHOOK) CALL DR_HOOK('FEMEANWS',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE FEMEANWS_WWM

! ----------------------------------------------------------------------

