      SUBROUTINE TAUHF(ML)

! ----------------------------------------------------------------------

!**** *TAUHF* - COMPUTATION OF HIGH-FREQUENCY STRESS.

!     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90

!*    PURPOSE.
!     ---------

!       COMPUTE HIGH-FREQUENCY WAVE STRESS

!**   INTERFACE.
!     ----------

!       *CALL* *TAUHF(ML)*
!             *ML*  NUMBER OF FREQUENCIES.

!     METHOD.
!     -------

!       SEE REFERENCE FOR WAVE STRESS CALCULATION.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FR
!      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,ALPHA    ,XKAPPA
!      USE YOWPCONS , ONLY : G        ,ZPI
!      USE YOWTABL  , ONLY : IUSTAR   ,IALPHA   ,USTARM   ,TAUHFT   ,
!     &            DELUST   ,DELALP

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, &
     &                      IUSTAR, IALPHA, USTARM, TAUHFT, STAT, &
     &                      DELUST, DELALP, ALPHA, BETAMAX, RKIND, &
     &                      XKAPPA, ZALP, &
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP, &
     &                      ZERO, ONE, &
     &                      SRCDBG

      IMPLICIT NONE

! ----------------------------------------------------------------------

      INTEGER, INTENT(IN)  :: ML
      INTEGER              :: I, J, K, L, M
      INTEGER, PARAMETER :: JTOT=250

      REAL(rkind), ALLOCATABLE :: W(:)
      REAL(rkind) :: ALPHAM, ALPHAMCOEF, CONST1, OMEGAC, X0, UST, Z0, OMEGACC, YC
      REAL(rkind) :: DELY, OMEGA, CM, ZX, ZARG, ZMU, ZLOG, Y, ZBETA
      integer istat

! ----------------------------------------------------------------------

!*    1. PRELIMINARY CALCULATIONS.
!        -------------------------

      ALPHAMCOEF = 40.
      ALPHAM = ALPHAMCOEF*ALPHA
      DELUST = USTARM/REAL(IUSTAR)
      DELALP = ALPHAM/REAL(IALPHA)

      CONST1 = BETAMAX/XKAPPA**2

      WRITE(100005,*) ALPHAMCOEF, ALPHAM, DELUST, DELALP, CONST1

      ALLOCATE(W(JTOT))
      W=1.
      W(1)=0.5
      W(JTOT)=0.5

      IF(.NOT.ALLOCATED(TAUHFT)) ALLOCATE(TAUHFT(0:IUSTAR,0:IALPHA,ML))

      DO M=1,ML

        OMEGAC = ZPI*FR(M)

        DO L=0,IALPHA
          DO K=0,IUSTAR
            TAUHFT(K,L,M) = 0.
          ENDDO
        ENDDO

!*    2. CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
!        ------------------------------------------------

        X0 = 0.05
        DO L=0,IALPHA
          DO K=0,IUSTAR
            UST      = MAX(REAL(K)*DELUST,0.000001)
            Z0       = UST**2*(ALPHA+REAL(L)*DELALP)/G
            OMEGACC  = MAX(OMEGAC,X0*G/UST)
            YC       = OMEGACC*SQRT(Z0/G)
            DELY     = MAX((1.-YC)/REAL(JTOT),0.)
            DO J=1,JTOT
              Y        = YC+REAL(J-1)*DELY
              OMEGA    = Y*SQRT(G/Z0)
              CM       = G/OMEGA
              ZX       = UST/CM +ZALP
              ZARG     = MIN(XKAPPA/ZX,20.)
              ZMU      = MIN(G*Z0/CM**2*EXP(ZARG),1.)

              ZLOG         = MIN(LOG(ZMU),0.)
              ZBETA        = CONST1*ZMU*ZLOG**4
              TAUHFT(K,L,M)= TAUHFT(K,L,M)+W(J)*ZBETA/Y*DELY
              !WRITE(100005,'(3I10,10F15.7)') L, K, J, DELY, ZMU, TAUHFT(K,L,M)
            ENDDO
          ENDDO
        ENDDO

      ENDDO

      DO M=1,ML
        WRITE(100000,*) M, FR(M)
      ENDDO

      DO M=1,ML
        DO L=0,IALPHA
          DO K=0,IUSTAR
            WRITE(100002,*) K,L,M,TAUHFT(K,L,M)
          ENDDO
        ENDDO
      ENDDO

      WRITE(111111,'(A10,F20.10)') 'TAUHFT', SUM(TAUHFT)

      DEALLOCATE(W)

      RETURN
      END SUBROUTINE TAUHF
     SUBROUTINE TAUHF_ECMWF_NEW

! ----------------------------------------------------------------------

!**** *TAUHF* - COMPUTATION OF HIGH-FREQUENCY STRESS.

!     PETER A.E.M. JANSSEN    KNMI      OCTOBER 90

!*    PURPOSE.
!     ---------

!       COMPUTE HIGH-FREQUENCY WAVE STRESS

!**   INTERFACE.
!     ----------

!       *CALL* *TAUHF(ML)*
!             *ML*  NUMBER OF FREQUENCIES.

!     METHOD.
!     -------

!       SEE REFERENCE FOR WAVE STRESS CALCULATION.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.

! ----------------------------------------------------------------------

!      USE YOWFRED  , ONLY : FR
!      USE YOWCOUP  , ONLY : BETAMAX  ,ZALP     ,ALPHA    ,XKAPPA
!      USE YOWPCONS , ONLY : G        ,ZPI
!      USE YOWTABL  , ONLY : IUSTAR   ,IALPHA   ,USTARM   ,TAUHFT   ,
!     &            DELUST   ,DELALP

       USE DATAPOOL, ONLY : FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, &
     &                      IUSTAR, IALPHA, USTARM, TAUHFT, STAT, &
     &                      DELUST, DELALP, ALPHA, BETAMAX, RKIND, &
     &                      XKAPPA, ZALP, &
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP, &
     &                      ZERO, ONE, & 
     &                      SRCDBG 

      IMPLICIT NONE

! ----------------------------------------------------------------------

      INTEGER              :: I, J, K, L, M
      INTEGER, PARAMETER :: JTOT=250

      REAL(rkind), ALLOCATABLE :: W(:)
      REAL(rkind) :: ALPHAM, ALPHAMCOEF, CONST1, OMEGAC, X0, UST, Z0, OMEGACC, YC
      REAL(rkind) :: DELY, OMEGA, CM, ZX, ZARG, ZMU, ZLOG, Y, ZBETA
      integer istat

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------

!*    1. PRELIMINARY CALCULATIONS.
!        -------------------------

      ALPHAMCOEF = 40.
      ALPHAM = ALPHAMCOEF*ALPHA
      DELUST = USTARM/REAL(IUSTAR)
      DELALP = ALPHAM/REAL(IALPHA)

      CONST1 = BETAMAX/XKAPPA**2

      WRITE(5011) DELUST, DELALP

      ALLOCATE(W(JTOT), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_ecmwf, allocate error 1')

      W=1.
      W(1)=0.5
      W(JTOT)=0.5

      DO M=1,NFRE

!        WRITE(STAT%FHNDL,*) 'DONE WITH M = ', M, 'OF   ', NFRE 

        OMEGAC = ZPI*FR(M)

        DO L=0,IALPHA
          DO K=0,IUSTAR
            TAUHFT(K,L,M) = 0.
          ENDDO
        ENDDO

!*    2. CALCULATE HIGH-FREQUENCY CONTRIBUTION TO STRESS.
!        ------------------------------------------------

        X0 = 0.05
        DO L=0,IALPHA
          DO K=0,IUSTAR
            UST      = MAX(REAL(K)*DELUST,0.000001_rkind)
            Z0       = UST**2*(ALPHA+REAL(L)*DELALP)/G
            OMEGACC  = MAX(OMEGAC,X0*G/UST)
            YC       = OMEGACC*SQRT(Z0/G)
            DELY     = MAX((1.-YC)/REAL(JTOT),ZERO)
            DO J=1,JTOT
              Y        = YC+REAL(J-1)*DELY
              OMEGA    = Y*SQRT(G/Z0)
              CM       = G/OMEGA
              ZX       = UST/CM +ZALP
              ZARG     = MIN(XKAPPA/ZX,20._rkind)
              ZMU      = MIN(G*Z0/CM**2*EXP(ZARG),ONE)

              ZLOG         = MIN(LOG(ZMU),ZERO)
              ZBETA        = CONST1*ZMU*ZLOG**4
              TAUHFT(K,L,M)= TAUHFT(K,L,M)+W(J)*ZBETA/Y*DELY
            ENDDO
          ENDDO
        ENDDO

      ENDDO

      WRITE(5011) TAUHFT

      DO M=1,NFRE
        WRITE(100000,*) M, FR(M)
      ENDDO

      DO M=1,NFRE
        DO L=0,IALPHA
          DO K=0,IUSTAR
            WRITE(100002,*) K,L,M,TAUHFT(K,L,M)
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(W)
      WRITE(SRCDBG%FHNDL,*) 'RESULTS HIGH FREQ STRESS TABLE'
      WRITE(SRCDBG%FHNDL,*) TAUHFT

      END SUBROUTINE TAUHF_ECMWF_NEW
