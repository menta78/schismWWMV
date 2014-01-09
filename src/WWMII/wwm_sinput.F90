      SUBROUTINE SINPUT (F, FL, IJS, IJL, THWNEW, USNEW, Z0NEW, &
     &                   ROAIRN, ZIDLNEW, SL, XLLWS)
! ----------------------------------------------------------------------

!**** *SINPUT* - COMPUTATION OF INPUT SOURCE FUNCTION.

!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990

!     OPTIMIZED BY : H. GUENTHER

!     MODIFIED BY : 
!       J-R BIDLOT NOVEMBER 1995
!       J-R BIDLOT FEBRUARY 1996-97
!       J-R BIDLOT FEBRUARY 1999 : INTRODUCE ICALL AND NCALL
!       P.A.E.M. JANSSEN MAY 2000 : INTRODUCE GUSTINESS
!       J-R BIDLOT FEBRUARY 2001 : MAKE IT FULLY IMPLICIT BY ONLY
!                                  USING NEW STRESS AND ROUGHNESS. 
!       S. ABDALLA OCTOBER 2001:  INTRODUCTION OF VARIABLE AIR
!                                 DENSITY AND STABILITY-DEPENDENT 
!                                 WIND GUSTINESS
!       P.A.E.M. JANSSEN OCTOBER 2008: INTRODUCE DAMPING WHEN WAVES ARE 
!                                      RUNNING FASTER THAN THE WIND.

!*    PURPOSE.
!     ---------

!       COMPUTE INPUT SOURCE FUNCTION AND STORE ADDITIVELY INTO NET
!       SOURCE FUNCTION ARRAY, ALSO COMPUTE FUNCTIONAL DERIVATIVE OF
!       INPUT SOURCE FUNCTION.
!
!       GUSTINESS IS INTRODUCED FOLL0WING THE APPROACH OF JANSSEN(1986),
!       USING A GAUSS-HERMITE APPROXIMATION SUGGESTED BY MILES(1997).
!       IN THE PRESENT VERSION ONLY TWO HERMITE POLYNOMIALS ARE UTILISED
!       IN THE EVALUATION OF THE PROBABILITY INTEGRAL. EXPLICITELY ONE THEN
!       FINDS:
!
!             <GAMMA(X)> = 0.5*( GAMMA(X(1+SIG)) + GAMMA(X(1-SIG)) )
!
!       WHERE X IS THE FRICTION VELOCITY AND SIG IS THE RELATIVE GUSTINESS
!       LEVEL.

!**   INTERFACE.
!     ----------

!     *CALL* *SINPUT (F, FL, IJS, IJL, THWNEW, USNEW, Z0NEW,
!    &                   ROAIRN,ZIDLNEW, SL, LLWS)
!            *F* - SPECTRUM.
!           *FL* - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!          *IJS* - INDEX OF FIRST GRIDPOINT.
!          *IJL* - INDEX OF LAST GRIDPOINT.
!       *THWNEW* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!        *USNEW* - NEW FRICTION VELOCITY IN M/S.
!        *Z0NEW* - ROUGHNESS LENGTH IN M.
!       *ROAIRN* - AIR DENSITY IN KG/M3
!      *ZIDLNEW* - Zi/L  USED FOR GUSTINESS.
!                  (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!           *SL* - TOTAL SOURCE FUNCTION ARRAY.
!         *LLWS* - TRUE WHERE SINPUT IS POSITIVE


!     METHOD.
!     -------

!       SEE REFERENCE.

!     EXTERNALS.
!     ----------

!       NONE.

!     MODIFICATIONS
!     -------------

!     - REMOVAL OF CALL TO CRAY SPECIFIC FUNCTIONS EXPHF AND ALOGHF
!       BY THEIR STANDARD FORTRAN EQUIVALENT EXP and ALOGHF
!     - MODIFIED TO MAKE INTEGRATION SCHEME FULLY IMPLICIT
!     - INTRODUCTION OF VARIABLE AIR DENSITY
!     - INTRODUCTION OF WIND GUSTINESS

!     REFERENCE.
!     ----------

!       P. JANSSEN, J.P.O., 1989.
!       P. JANSSEN, J.P.O., 1991

! ----------------------------------------------------------------------

      !USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,XKAPPA
      !USE YOWFRED  , ONLY : FR       ,TH
      !USE YOWPARAM , ONLY : NANG     ,NFRE
      !USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER   ,YEPS
      !USE YOWSHAL  , ONLY : TFAK     ,INDEP
      !USE YOWSTAT  , ONLY : ISHALLO  ,IDAMPING
      !USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, RKIND, &
     &                DFIM, DFIMOFR, DFFR, DFFR2, WK, ZALP, &
     &                IUSTAR, IALPHA, USTARM, TAUT, ONETHIRD, RKIND, &
     &                DELUST, DELALP, BETAMAX, XKAPPA, IDAMPING, &
     &                ROWATER => RHOW, &
     &                TH => SPDIR, &
     &                DELTH => DDIR, &
     &                G => G9, &
     &                ZPI => PI2, &
     &                EPSMIN => SMALL, &
     &                NANG => MDC, &
     &                NFRE => MSC, &
     &                INDEP => DEP
      IMPLICIT NONE

! ----------------------------------------------------------------------

!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      REAL(rkind),DIMENSION(IJS:IJL,NANG,NFRE) :: F, FL, SL
      REAL(rkind),DIMENSION(IJS:IJL) :: THWNEW, USNEW, Z0NEW, ROAIRN, ZIDLNEW

! ----------------------------------------------------------------------

      REAL(rkind), DIMENSION(IJS:IJL,NANG) :: TEMP1, UFAC2
      REAL(rkind), DIMENSION(IJS:IJL) :: UCN1, UCN2, ZCN, CM, USP, USM
      REAL(rkind), DIMENSION(IJS:IJL) :: SIG_N, XV1D, XV2D
      REAL(rkind), DIMENSION(IJS:IJL) :: CNSN
      REAL(rkind), DIMENSION(NFRE) :: FAC, CONST
      REAL(rkind), DIMENSION(IJS:IJL,NANG,NFRE) :: XLLWS
      LOGICAL L1,L2,LZ(IJS:IJL,NANG)
      REAL(rkind) X1,X2,X1D,X2D,ZLOG1,ZLOG2,CONST3,XV1,XV2,ZBETA1,ZBETA2
      REAL(rkind) XKAPPAD, BG_GUST, CONST1, DC_DDU, C_D, U10, SIG_CONV
      REAL(rkind) TEMPD(IJS:IJL,NANG),UCN1D(IJS:IJL),UCN2D(IJS:IJL), ZLOG2X
      REAL, PARAMETER :: A = 0.8d0/1000.d0
      REAL, PARAMETER :: B = 0.08d0/1000.d0
      INTEGER :: I, J, K, IJ, M, IJS, IJL
      !REAL ZHOOK_HANDLE

      !IF (LHOOK) CALL DR_HOOK('SINPUT',0,ZHOOK_HANDLE)

! ----------------------------------------------------------------------

      BG_GUST  = 0.d0        ! NO BACKGROUND GUSTINESS (S0 12. IS NOT USED)
      CONST1   = BETAMAX/XKAPPA**2 
      CONST3   = 2.*XKAPPA/CONST1  ! SEE IDAMPING
      XKAPPAD  = 1.D0/XKAPPA


      CONST3 = IDAMPING*CONST3

      WRITE(111116, '(10F15.7)') CONST1, CONST3, XKAPPAD, CONST3 

!*    1. PRECALCULATED ANGULAR DEPENDENCE.
!        ---------------------------------

      DO K=1,NANG
        DO IJ=IJS,IJL
          TEMP1(IJ,K) = COS(TH(K)-THWNEW(IJ))
          IF(TEMP1(IJ,K) .GT. 0.01) THEN
            LZ(IJ,K) = .TRUE.
            TEMPD(IJ,K) = 1.D0/TEMP1(IJ,K)
          ELSE
            LZ(IJ,K) = .FALSE.
            TEMPD(IJ,K) = 1.D0
          ENDIF
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
!
!       IN THE FOLLOWING U10 IS ESTIMATED ASSUMING EVERYTHING IS
!       BASED ON U*
!
        U10 = USNEW(IJ)/XKAPPA*LOG(10./Z0NEW(IJ))
        C_D = A+B*U10
        DC_DDU = B
        SIG_CONV = 1. + 0.5*U10/C_D*DC_DDU
        SIG_N (IJ) = MIN(0.5, SIG_CONV * &
     &                        (BG_GUST*USNEW(IJ)**3+ &
     &                         0.5*XKAPPA*ZIDLNEW(IJ)**3)**ONETHIRD &
     &                           /U10 &
     &                 )
        WRITE(111116, '(I10,10F15.7)') IJ, U10, C_D, DC_DDU, SIG_CONV, SIG_N(IJ), USP(IJ), USM(IJ)
      ENDDO

      DO IJ=IJS,IJL
        USP(IJ) = USNEW(IJ)*(1.+SIG_N(IJ))
        USM(IJ) = USNEW(IJ)*(1.-SIG_N(IJ))
      ENDDO

! ----------------------------------------------------------------------

!*    2. LOOP OVER FREQUENCIES.
!        ----------------------


      DO M=1,NFRE

        FAC(M) = ZPI*FR(M)
        CONST(M)=FAC(M)*CONST1

!*      INVERSE OF PHASE VELOCITIES.
!       ----------------------------

        IF (ISHALLO.EQ.1) THEN
          DO IJ=IJS,IJL
            CM(IJ) = FAC(M)/G
          ENDDO
        ELSE
          DO IJ=IJS,IJL
            !CM(IJ) = TFAK(INDEP(IJ),M)/FAC(M)
            CM(IJ) = WK(IJ,M)/FAC(M)
          ENDDO
        ENDIF

!*      PRECALCULATE FREQUENCY DEPENDENCE.
!       ----------------------------------

        DO IJ=IJS,IJL
          UCN1(IJ) = USP(IJ)*CM(IJ) + ZALP
          UCN2(IJ) = USM(IJ)*CM(IJ) + ZALP

          UCN1D(IJ) = 1.D0/ UCN1(IJ)
          UCN2D(IJ) = 1.D0/ UCN2(IJ)

          ZCN(IJ)  = LOG(G*Z0NEW(IJ)*CM(IJ)**2)
          CNSN(IJ) = CONST(M) * ROAIRN(IJ)/ROWATER

          XV1      = -USP(IJ)*XKAPPAD*ZCN(IJ)*CM(IJ)
          XV2      = -USM(IJ)*XKAPPAD*ZCN(IJ)*CM(IJ)

          XV1D(IJ) = 1.D0/XV1
          XV2D(IJ) = 1.D0/XV2

        ENDDO

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

        DO K=1,NANG
          DO IJ=IJS,IJL
            X1    = TEMP1(IJ,K)*UCN1(IJ)
            X1D   = TEMPD(IJ,K)*UCN1D(IJ)
            ZLOG1 = ZCN(IJ) + XKAPPA*X1D
            L1    = ZLOG1.LT.0.
            X2    = TEMP1(IJ,K)*UCN2(IJ)
            X2D   = TEMPD(IJ,K)*UCN2D(IJ)
            ZLOG2 = ZCN(IJ) + XKAPPA*X2D
            L2    = ZLOG2.LT.0.

            ZBETA1 = CONST3*(TEMP1(IJ,K)-XV1D(IJ))*UCN1(IJ)**2 
            ZBETA2 = CONST3*(TEMP1(IJ,K)-XV2D(IJ))*UCN2(IJ)**2             

            IF (LZ(IJ,K)) THEN
              IF (L1) THEN
                ZLOG2X=ZLOG1*ZLOG1*X1
                UFAC2(IJ,K) = EXP(ZLOG1)*ZLOG2X*ZLOG2X+ZBETA1
                XLLWS(IJ,K,M)= 1.
              ELSE
                UFAC2(IJ,K) = ZBETA1
                XLLWS(IJ,K,M)= 0.
              ENDIF
              IF (L2) THEN
                ZLOG2X=ZLOG2*ZLOG2*X2
                UFAC2(IJ,K) = UFAC2(IJ,K)+  &
     &                        EXP(ZLOG2)*ZLOG2X*ZLOG2X+ZBETA2
                XLLWS(IJ,K,M)= 1.
              ELSE
                UFAC2(IJ,K) = UFAC2(IJ,K)+ZBETA2
              ENDIF
            ELSE
              UFAC2(IJ,K) = ZBETA1+ZBETA2
              XLLWS(IJ,K,M)= 0.
            ENDIF
            WRITE(111116, '(3I10,10F15.7)') M, K, IJ, UFAC2(IJ,M), XLLWS(IJ, K,M), FAC(M), CONST(M), CM(IJ), XV1D(IJ), XV2D(IJ)
          ENDDO
        ENDDO

!*    2.2 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
!         ------------------------------------------------

        DO K=1,NANG
          DO IJ=IJS,IJL
            FL(IJ,K,M) = 0.5*CNSN(IJ)*UFAC2(IJ,K)
            SL(IJ,K,M) = FL(IJ,K,M)*F(IJ,K,M)
      WRITE(111116, '(3I10,10F15.7)') M, K, IJ, FL(IJ,K,M), SL(IJ,K,M)
          ENDDO
        ENDDO

      ENDDO


      !IF (LHOOK) CALL DR_HOOK('SINPUT',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE SINPUT

      SUBROUTINE SINPUT_WWM (IP, F, FL, THWNEW, USNEW, Z0NEW, ROAIRN, ZIDLNEW, SL, XLLWS)
! ----------------------------------------------------------------------

!**** *SINPUT* - COMPUTATION OF INPUT SOURCE FUNCTION.

!     P.A.E.M. JANSSEN    KNMI      AUGUST    1990

!     OPTIMIZED BY : H. GUENTHER

!     MODIFIED BY : 
!       J-R BIDLOT NOVEMBER 1995
!       J-R BIDLOT FEBRUARY 1996-97
!       J-R BIDLOT FEBRUARY 1999 : INTRODUCE ICALL AND NCALL
!       P.A.E.M. JANSSEN MAY 2000 : INTRODUCE GUSTINESS
!       J-R BIDLOT FEBRUARY 2001 : MAKE IT FULLY IMPLICIT BY ONLY
!                                  USING NEW STRESS AND ROUGHNESS. 
!       S. ABDALLA OCTOBER 2001:  INTRODUCTION OF VARIABLE AIR
!                                 DENSITY AND STABILITY-DEPENDENT 
!                                 WIND GUSTINESS
!       P.A.E.M. JANSSEN OCTOBER 2008: INTRODUCE DAMPING WHEN WAVES ARE 
!                                      RUNNING FASTER THAN THE WIND.

!*    PURPOSE.
!     ---------

!       COMPUTE INPUT SOURCE FUNCTION AND STORE ADDITIVELY INTO NET
!       SOURCE FUNCTION ARRAY, ALSO COMPUTE FUNCTIONAL DERIVATIVE OF
!       INPUT SOURCE FUNCTION.
!
!       GUSTINESS IS INTRODUCED FOLL0WING THE APPROACH OF JANSSEN(1986),
!       USING A GAUSS-HERMITE APPROXIMATION SUGGESTED BY MILES(1997).
!       IN THE PRESENT VERSION ONLY TWO HERMITE POLYNOMIALS ARE UTILISED
!       IN THE EVALUATION OF THE PROBABILITY INTEGRAL. EXPLICITELY ONE THEN
!       FINDS:
!
!             <GAMMA(X)> = 0.5*( GAMMA(X(1+SIG)) + GAMMA(X(1-SIG)) )
!
!       WHERE X IS THE FRICTION VELOCITY AND SIG IS THE RELATIVE GUSTINESS
!       LEVEL.

!**   INTERFACE.
!     ----------

!     *CALL* *SINPUT (F, FL, IJS, IJL, THWNEW, USNEW, Z0NEW,
!    &                   ROAIRN,ZIDLNEW, SL, LLWS)
!            *F* - SPECTRUM.
!           *FL* - DIAGONAL MATRIX OF FUNCTIONAL DERIVATIVE.
!          *IJS* - INDEX OF FIRST GRIDPOINT.
!          *IJL* - INDEX OF LAST GRIDPOINT.
!       *THWNEW* - WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                  NOTATION (POINTING ANGLE OF WIND VECTOR,
!                  CLOCKWISE FROM NORTH).
!        *USNEW* - NEW FRICTION VELOCITY IN M/S.
!        *Z0NEW* - ROUGHNESS LENGTH IN M.
!       *ROAIRN* - AIR DENSITY IN KG/M3
!      *ZIDLNEW* - Zi/L  USED FOR GUSTINESS.
!                  (Zi: INVERSION HEIGHT, L: MONIN-OBUKHOV LENGTH).
!           *SL* - TOTAL SOURCE FUNCTION ARRAY.
!         *LLWS* - TRUE WHERE SINPUT IS POSITIVE


!     METHOD.
!     -------

!       SEE REFERENCE.

!     EXTERNALS.
!     ----------

!       NONE.

!     MODIFICATIONS
!     -------------

!     - REMOVAL OF CALL TO CRAY SPECIFIC FUNCTIONS EXPHF AND ALOGHF
!       BY THEIR STANDARD FORTRAN EQUIVALENT EXP and ALOGHF
!     - MODIFIED TO MAKE INTEGRATION SCHEME FULLY IMPLICIT
!     - INTRODUCTION OF VARIABLE AIR DENSITY
!     - INTRODUCTION OF WIND GUSTINESS

!     REFERENCE.
!     ----------

!       P. JANSSEN, J.P.O., 1989.
!       P. JANSSEN, J.P.O., 1991

! ----------------------------------------------------------------------

      !USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,XKAPPA
      !USE YOWFRED  , ONLY : FR       ,TH
      !USE YOWPARAM , ONLY : NANG     ,NFRE
      !USE YOWPCONS , ONLY : G        ,ZPI      ,ROWATER   ,YEPS
      !USE YOWSHAL  , ONLY : TFAK     ,INDEP
      !USE YOWSTAT  , ONLY : ISHALLO  ,IDAMPING
      !USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      !USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPRM

      USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, &
     &                DFIM, DFIMOFR, DFFR, DFFR2, WK, ZALP, &
     &                IUSTAR, IALPHA, USTARM, TAUT, ONETHIRD, RKIND, &
     &                DELUST, DELALP, BETAMAX, XKAPPA, IDAMPING, &
     &                ROWATER => RHOW, &
     &                TH => SPDIR, &
     &                DELTH => DDIR, &
     &                G => G9, &
     &                ZPI => PI2, &
     &                EPSMIN => SMALL, &
     &                NANG => MDC, &
     &                NFRE => MSC, &
     &                INDEP => DEP
      IMPLICIT NONE
! ----------------------------------------------------------------------

!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      INTEGER, INTENT(IN) :: IP

      REAL(rkind),DIMENSION(NANG,NFRE) :: F, FL, SL
      REAL(rkind) :: THWNEW, USNEW, Z0NEW, ROAIRN, ZIDLNEW

! ----------------------------------------------------------------------

      REAL(rkind), DIMENSION(NANG) :: TEMP1, UFAC2
      REAL(rkind) :: UCN1, UCN2, ZCN, CM, USP, USM
      REAL(rkind) :: SIG_N, XV1D, XV2D
      REAL(rkind) :: CNSN
      REAL(rkind), DIMENSION(NFRE) :: FAC, CONST
      REAL(rkind), DIMENSION(NANG,NFRE) :: XLLWS
      LOGICAL :: L1,L2,LZ(NANG)
      REAL(rkind) :: X1,X2,X1D,X2D,ZLOG1,ZLOG2,CONST3,XV1,XV2,ZBETA1,ZBETA2
      REAL(rkind) :: XKAPPAD, BG_GUST, CONST1, U10, C_D, DC_DDU, SIG_CONV
      REAL(rkind) :: ZLOG2X
      INTEGER :: IJ, K, L, M, N 
      REAL(rkind) TEMPD(NANG),UCN1D,UCN2D
      REAL(rkind), PARAMETER :: A = 0.8/1000. 
      REAL(rkind), PARAMETER :: B = 0.08/1000.
      !REAL(KIND=JPRB) ZHOOK_HANDLE

      !IF (LHOOK) CALL DR_HOOK('SINPUT',0,ZHOOK_HANDLE)

! ----------------------------------------------------------------------

      BG_GUST  = 0.        ! NO BACKGROUND GUSTINESS (S0 12. IS NOT USED)
      CONST1   = BETAMAX/XKAPPA**2
      CONST3   = 2.*XKAPPA/CONST1  ! SEE IDAMPING
      XKAPPAD  = 1.D0/XKAPPA

      CONST3 = IDAMPING*CONST3

!*    1. PRECALCULATED ANGULAR DEPENDENCE.
!        ---------------------------------

      DO K=1,NANG
        TEMP1(K) = COS(TH(K)-THWNEW)
        IF(TEMP1(K) .GT. 0.01) THEN
          LZ(K) = .TRUE.
          TEMPD(K) = 1.D0/TEMP1(K)
        ELSE
          LZ(K) = .FALSE.
          TEMPD(K) = 1.D0
        ENDIF
      ENDDO

!
!       IN THE FOLLOWING U10 IS ESTIMATED ASSUMING EVERYTHING IS
!       BASED ON U*
!
      U10 = USNEW/XKAPPA*LOG(10./Z0NEW)
      C_D = A+B*U10
      DC_DDU = B
      SIG_CONV = 1. + 0.5*U10/C_D*DC_DDU
      SIG_N = MIN(0.5_rkind, SIG_CONV * (BG_GUST*USNEW**3+ 0.5*XKAPPA*ZIDLNEW**3)**ONETHIRD/U10)


      USP = USNEW*(1.+SIG_N)
      USM = USNEW*(1.-SIG_N)

      WRITE(111116, '(10F15.7)') U10, C_D, DC_DDU, SIG_CONV, SIG_N, USP, USM
! ----------------------------------------------------------------------

!*    2. LOOP OVER FREQUENCIES.
!        ----------------------


      DO M=1,NFRE

        FAC(M) = ZPI*FR(M)
        CONST(M)=FAC(M)*CONST1

!*      INVERSE OF PHASE VELOCITIES.
!       ----------------------------

        IF (ISHALLO.EQ.1) THEN
          CM = FAC(M)/G
        ELSE
          !CM(IJ) = TFAK(INDEP(IJ),M)/FAC(M)
          CM = WK(IP,M)/FAC(M)
        ENDIF

!*      PRECALCULATE FREQUENCY DEPENDENCE.
!       ----------------------------------
        UCN1 = USP*CM + ZALP
        UCN2 = USM*CM + ZALP
        UCN1D = 1.D0/ UCN1
        UCN2D = 1.D0/ UCN2
        ZCN  = LOG(G*Z0NEW*CM**2)
        CNSN = CONST(M) * ROAIRN/ROWATER
        XV1      = -USP*XKAPPAD*ZCN*CM
        XV2      = -USM*XKAPPAD*ZCN*CM
        XV1D = 1.D0/XV1
        XV2D = 1.D0/XV2

        WRITE(111116, '(I10,10F15.7)') M, FAC(M), CONST(M), CM, XV1D, XV2D 

!*    2.1 LOOP OVER DIRECTIONS.
!         ---------------------

        DO K=1,NANG
          X1    = TEMP1(K)*UCN1
          X1D   = TEMPD(K)*UCN1D
          ZLOG1 = ZCN + XKAPPA*X1D
          L1    = ZLOG1.LT.0.
          X2    = TEMP1(K)*UCN2
          X2D   = TEMPD(K)*UCN2D
          ZLOG2 = ZCN + XKAPPA*X2D
          L2    = ZLOG2.LT.0.

          ZBETA1 = CONST3*(TEMP1(K)-XV1D)*UCN1**2
          ZBETA2 = CONST3*(TEMP1(K)-XV2D)*UCN2**2
          IF (LZ(K)) THEN
            IF (L1) THEN
              ZLOG2X=ZLOG1*ZLOG1*X1
              UFAC2(K) = EXP(ZLOG1)*ZLOG2X*ZLOG2X+ZBETA1
              XLLWS(K,M)= 1.
            ELSE
              UFAC2(K) = ZBETA1
              XLLWS(K,M)= 0.
            ENDIF
            IF (L2) THEN
              ZLOG2X=ZLOG2*ZLOG2*X2
              UFAC2(K) = UFAC2(K)+EXP(ZLOG2)*ZLOG2X*ZLOG2X+ZBETA2
              XLLWS(K,M)= 1.
            ELSE
              UFAC2(K) = UFAC2(K)+ZBETA2
            ENDIF
          ELSE
            UFAC2(K) = ZBETA1+ZBETA2
            XLLWS(K,M)= 0.
          ENDIF

          WRITE(111116, '(I10,10F15.7)') M, K, L1, L2, UFAC2(M), ZBETA2, XLLWS(K,M)

        ENDDO

!*    2.2 ADDING INPUT SOURCE TERM TO NET SOURCE FUNCTION.
!         ------------------------------------------------

        DO K=1,NANG
          FL(K,M) = 0.5*CNSN*UFAC2(K)
          SL(K,M) = FL(K,M)*F(K,M)
          !write(DBG%FHNDL,'(2I10,4F15.8)') M, K, FL(K,M), SL(K,M), UFAC2(K), F(K,M)
        ENDDO

      ENDDO


!      IF (LHOOK) CALL DR_HOOK('SINPUT',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE SINPUT_WWM

