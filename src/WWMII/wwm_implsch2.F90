! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
    SUBROUTINE PREINTRHS (FL3, FL, IJS, IJL, IG, &
     &                    THWOLD, USOLD, &
     &                    TAUW, Z0OLD, &
     &                    ROAIRO, ZIDLOLD, &
     &                    U10NEW, THWNEW, USNEW, &
     &                    Z0NEW, ROAIRN, ZIDLNEW, &
     &                    SL, FCONST, FMEANWS, MIJ)

       USE DATAPOOL, ONLY : MNP, FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, FRINTF, COFRM4, CG, WK, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, EMEAN, FMEAN, TH, RKIND, DELU, &
     &                      JUMAX, DT4S, FRM5, IPHYS, LOUTWAM, CD, UFRIC, ALPHA_CH, Z0, ITEST, LCFLX, &
     &                      DELTH => DDIR, TESTNODE, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP, &
     &                      ROWATER => RHOW, &
     &                      RHOAIR => RHOA

      IMPLICIT NONE

! ----------------------------------------------------------------------
!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL(rkind) :: FL(IJS:IJL,NANG,NFRE),FL3(IJS:IJL,NANG,NFRE),SL(IJS:IJL,NANG,NFRE)
      REAL(rkind),DIMENSION(IJS:IJL,NFRE) :: FCONST
      REAL(rkind),DIMENSION(IJS:IJL) :: THWOLD,USOLD,Z0OLD,TAUW, &
     &                           ROAIRO,ZIDLOLD,FMEANWS 
      REAL(rkind),DIMENSION(IJS:IJL) :: U10NEW,THWNEW,USNEW,Z0NEW, &
     &                           ROAIRN,ZIDLNEW

! ----------------------------------------------------------------------
 
      INTEGER :: IJ,IJS,IJL,K,L,M,IG,ILEV,IDELT,IU06
      INTEGER :: JU(IJS:IJL)
      INTEGER :: MIJ(IJS:IJL)
      REAL(rkind) :: GTEMP1, GTEMP2, FLHAB, XJ, DELT, DELT5, XIMP, AKM1
      REAL(rkind) :: AK2VGM1, XN, PHIDIAG, TAU
      REAL(rkind), DIMENSION(NFRE) :: DELFL
      REAL(rkind), DIMENSION(IJS:IJL) :: EMEANWS, USFM, GADIAG 
      REAL(rkind), DIMENSION(IJS:IJL) :: F1MEAN, AKMEAN, XKMEAN
      REAL(rkind), DIMENSION(IJS:IJL) :: PHIEPS, TAUOC, PHIAW
      REAL(rkind), DIMENSION(IJS:IJL) :: TAUWLF,TAUWD,PHIAWDIAG,PHIAWUNR,PHIOC,PHIWA 
      REAL(rkind), DIMENSION(IJS:IJL,NANG) :: SPRD
      REAL(rkind), DIMENSION(IJS:IJL,NFRE) :: TEMP, TEMP2
      REAL(rkind), DIMENSION(IJS:IJL,NANG,NFRE) :: CIWAB 
      REAL(rkind), DIMENSION(IJS:IJL,NANG,NFRE) :: XLLWS

      IDELT = INT(DT4S)

!      REAL ZHOOK_HANDLE

!      IF (LHOOK) CALL DR_HOOK('IMPLSCH',0,ZHOOK_HANDLE)
! ----------------------------------------------------------------------

!*    1. INITIALISATION.
!        ---------------

      ! LCFLX=LWFLUX.OR.LWFLUXOUT.OR.LWNEMOCOU
! ----------------------------------------------------------------------

!*    2. COMPUTATION OF IMPLICIT INTEGRATION.
!        ------------------------------------

!         INTEGRATION IS DONE FROM CDATE UNTIL CDTPRO FOR A BLOCK
!         OF LATITUDES BETWEEN PROPAGATION CALLS.


!*    2.2 COMPUTE MEAN PARAMETERS.
!        ------------------------

      CALL FKMEAN(FL3, IJS, IJL, EMEAN(IJS), FMEAN(IJS), &
     &            F1MEAN, AKMEAN, XKMEAN)

      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'HS and TM'
      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(10F15.7)') 4*SQRT(EMEAN(IJS)), FMEAN(IJS), SUM(FL3)
      !WRITE(55555,*) 4*SQRT(EMEAN(IJS)), FMEAN(IJS)

      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'DIRECTIONAL PROPERTIES'
      DO K=1,NANG
        DO IJ=IJS,IJL
          SPRD(IJ,K)=MAX(0.,COS(TH(K)-THWNEW(IJ)))**2
!          WRITE(111113,'(I10,10F15.7)') K, SPRD(IJ,K), TH(K), THWNEW(IJ) 
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        XJ=U10NEW(IJ)/DELU
        JU(IJ)=MIN(JUMAX, MAX(NINT(XJ),1))
      ENDDO

! ----------------------------------------------------------------------

!*    2.3 COMPUTATION OF SOURCE FUNCTIONS.
!         --------------------------------

!*      2.3.1 INITIALISE SOURCE FUNCTION AND DERIVATIVE ARRAY.
!             ------------------------------------------------

!!            FL AND SL ARE INITIALISED IN SINPUT

      
      ILEV=1
      CALL AIRSEA (U10NEW(IJS), TAUW(IJS), USNEW(IJS), Z0NEW(IJS), &
     &   IJS, IJL, ILEV)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED BEFORE DO LOOP'
        CALL FLUSH (IU06)
      ENDIF

      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER AIRSEA 1'
      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,10F15.7)') IJS, U10NEW(IJS), TAUW(IJS), &
      &                              USNEW(IJS), Z0NEW(IJS), ILEV

!*    2.3.2 ADD SOURCE FUNCTIONS AND WAVE STRESS.
!           -------------------------------------

      IF(IPHYS.EQ.0) THEN
        CALL SINPUT (FL3, FL, IJS, IJL, THWNEW, USNEW, Z0NEW, &
     &             ROAIRN, ZIDLNEW, SL, XLLWS)
      ELSE
        CALL SINPUT_ARD (FL3, FL, IJS, IJL, THWNEW, USNEW, Z0NEW, &
     &             ROAIRN, ZIDLNEW, SL, XLLWS)
      ENDIF
      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER SINPUT 1'
      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(SL)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SINPUT CALLED'
        CALL FLUSH (IU06)
      ENDIF

!     MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
      CALL FEMEANWS(FL3,IJS,IJL,EMEANWS,FMEANWS,XLLWS)

!     COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
      CALL FRCUTINDEX(IJS, IJL, FMEAN(IJS), FMEANWS, MIJ)

      CALL STRESSO (FL3, IJS, IJL, THWNEW, USNEW, Z0NEW, &
     &              ROAIRN, TAUW, TAUWLF, PHIWA, &
     &              PHIAWDIAG, PHIAWUNR, SL, &
     &              MIJ, LCFLX)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: STRESSO CALLED'
        CALL FLUSH (IU06)
      ENDIF

      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER STRESSO 1'
      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(2I10,15F15.7)') IJS, IJL, SUM(FL3), &
     &              THWNEW, USNEW, Z0NEW, &
     &              ROAIRN, TAUW, TAUWLF, PHIWA, &
     &              PHIAWDIAG, PHIAWUNR, SUM(SL), &
     &              MIJ(IJS)

      CALL AIRSEA (U10NEW(IJS), TAUW(IJS), USNEW(IJS), Z0NEW(IJS), &
     &             IJS, IJL, ILEV)
      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER AIRSEA 2'
      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,10F15.7)') IJS, U10NEW(IJS), TAUW(IJS), &
     &             USNEW(IJS), Z0NEW(IJS)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED'
        CALL FLUSH (IU06)
      ENDIF

!     2.3.3 ADD THE OTHER SOURCE TERMS.
!           ---------------------------
      CALL SNONLIN (FL3, FL, IJS, IJL, IG, SL, AKMEAN(IJS))
      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER SNON'
      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(SL)
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SNONLIN CALLED'
        CALL FLUSH (IU06)
      ENDIF
      IF(IPHYS.EQ.0) THEN
        CALL SDISSIP (FL3 ,FL, IJS, IJL, IG, SL, F1MEAN, XKMEAN,&
     &                PHIOC, TAUWD, MIJ)
      ELSE
        CALL SDISS_ARDH_VEC (FL3 ,FL, IJS, IJL, SL, F1MEAN, XKMEAN,&
     &                PHIOC, TAUWD, MIJ)
      ENDIF
      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER DISSIP' 
      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(SL) 
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SDISSIP CALLED'
        CALL FLUSH (IU06)
      ENDIF
!SHALLOW
      IF(ISHALLO.NE.1) CALL SBOTTOM (FL3, FL, IJS, IJL, IG, SL)
!SHALLOW
      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER SBOTTOM' 
      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(SL) 

      DO IJ=IJS,IJL
        USOLD(IJ) = USNEW(IJ)
        Z0OLD(IJ) = Z0NEW(IJ)
        ROAIRO(IJ) = ROAIRN(IJ)
        ZIDLOLD(IJ) = ZIDLNEW(IJ)
      ENDDO

      RETURN
      END SUBROUTINE PREINTRHS 
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
    SUBROUTINE INTSPECWAM (FL3, FL, IJS, IJL, IG, &
     &                    THWOLD, USOLD, &
     &                    TAUW, Z0OLD, &
     &                    ROAIRO, ZIDLOLD, &
     &                    U10NEW, THWNEW, USNEW, &
     &                    Z0NEW, ROAIRN, ZIDLNEW, &
     &                    SL, FCONST, FMEANWS, MIJ)

       USE DATAPOOL, ONLY : MNP, FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, FRINTF, COFRM4, CG, WK, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, EMEAN, FMEAN, TH, RKIND, DELU, &
     &                      JUMAX, DT4S, FRM5, IPHYS, LOUTWAM, CD, UFRIC, ALPHA_CH, Z0, ITEST, LCFLX, &
     &                      DELTH => DDIR, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP, &
     &                      ROWATER => RHOW, &
     &                      RHOAIR => RHOA

      IMPLICIT NONE

! ----------------------------------------------------------------------
!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL(rkind) :: FL(IJS:IJL,NANG,NFRE),FL3(IJS:IJL,NANG,NFRE),SL(IJS:IJL,NANG,NFRE)
      REAL(rkind),DIMENSION(IJS:IJL,NFRE) :: FCONST
      REAL(rkind),DIMENSION(IJS:IJL) :: THWOLD,USOLD,Z0OLD,TAUW, &
     &                           ROAIRO,ZIDLOLD
      REAL(rkind),DIMENSION(IJS:IJL) :: U10NEW,THWNEW,USNEW,Z0NEW, &
     &                           ROAIRN,ZIDLNEW

! ----------------------------------------------------------------------

      INTEGER :: IJ,IJS,IJL,K,L,M,IG,ILEV,IDELT,IU06
      INTEGER :: JU(IJS:IJL)
      INTEGER :: MIJ(IJS:IJL)
      REAL(rkind) :: GTEMP1, GTEMP2, FLHAB, XJ, DELT, DELT5, XIMP, AKM1
      REAL(rkind) :: AK2VGM1, XN, PHIDIAG, TAU
      REAL(rkind), DIMENSION(NFRE) :: DELFL
      REAL(rkind), DIMENSION(IJS:IJL) :: EMEANWS, USFM, GADIAG
      REAL(rkind), DIMENSION(IJS:IJL) :: F1MEAN, AKMEAN, XKMEAN
      REAL(rkind), DIMENSION(IJS:IJL) :: PHIEPS, TAUOC, PHIAW, FMEANWS
      REAL(rkind), DIMENSION(IJS:IJL) :: TAUWLF,TAUWD,PHIAWDIAG,PHIAWUNR,PHIOC,PHIWA
      REAL(rkind), DIMENSION(IJS:IJL,NANG) :: SPRD
      REAL(rkind), DIMENSION(IJS:IJL,NFRE) :: TEMP, TEMP2
      REAL(rkind), DIMENSION(IJS:IJL,NANG,NFRE) :: CIWAB
      REAL(rkind), DIMENSION(IJS:IJL,NANG,NFRE) :: XLLWS

      IDELT = INT(DT4S)

!*    2.4 COMPUTATION OF NEW SPECTRA.
!         ---------------------------

!       INCREASE OF SPECTRUM IN A TIME STEP IS LIMITED TO A FINITE
!       FRACTION OF A TYPICAL F**(-4) EQUILIBRIUM SPECTRUM.

      DELT = IDELT
      XIMP = 1.0
      DELT5 = XIMP*DELT

      DO M=1,NFRE
        DELFL(M) = COFRM4(M)*DELT
      ENDDO
      DO IJ=IJS,IJL
        USFM(IJ) = USNEW(IJ)*MAX(FMEANWS(IJ),FMEAN(IJ))
        !IF (LOUTWAM) WRITE(111113,'(4F20.10)') USNEW(IJ), FMEANWS(IJ), FMEAN(IJ)
      ENDDO

      DO M=1,NFRE
        DO IJ=IJS,IJL
          TEMP(IJ,M) = USFM(IJ)*DELFL(M)
          !WRITE(111113,'(4F20.10)') DELFL(M), COFRM4(M), DELT
        ENDDO
      ENDDO

!      IF(ISHALLO.EQ.1) THEN
!        DO M=1,NFRE
!          DO IJ=IJS,IJL
!            TEMP2(IJ,M) = FRM5(M)
!          ENDDO
!        ENDDO
!      ELSE
!        DO M=1,NFRE
!          DO IJ=IJS,IJL
!AR: WAM TABLE REPLACES BY WWM WK            AKM1 = 1./TFAK(INDEP(IJ),M)
!            AKM1 = 1./WK(IJ,M)
!AR: WAM TABLE REPLACES BY WWM CG            AK2VGM1 = AKM1**2/TCGOND(INDEP(IJ),M)
!            AK2VGM1 = AKM1**2/CG(IJ,M)
!            TEMP2(IJ,M) = AKM1*AK2VGM1
!            WRITE(111113,'(4F20.10)') AKM1, AK2VGM1, TEMP2(IJ,M) 
!          ENDDO
!        ENDDO
!      ENDIF
!      WRITE(111113,*) 'MORE TEST'
      DO K=1,NANG
        DO M=1,NFRE
          DO IJ=IJS,IJL
            GTEMP1 = MAX((1.-DELT5*FL(IJ,K,M)),1.)
            GTEMP2 = DELT*SL(IJ,K,M)/GTEMP1
            FLHAB = ABS(GTEMP2)
            FLHAB = MIN(FLHAB,TEMP(IJ,M))
            FL3(IJ,K,M) = FL3(IJ,K,M) + SIGN(FLHAB,GTEMP2) 
!AR: ICE            FLLOWEST = FLMINFR(JU(IJ),M)*SPRD(IJ,K)
!AR: ICE            FL3(IJ,K,M) = MAX(FL3(IJ,K,M),FLLOWEST)
!      WRITE(111113,'(4F20.10)')GTEMP2,FLHAB,TEMP(IJ,M)
          ENDDO
        ENDDO
      ENDDO

      !IF (LOUTWAM) WRITE(111113,*) 'AFTER INTEGRATION' 
      !IF (LOUTWAM) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(SL), SUM(FL3)

      RETURN
      END SUBROUTINE INTSPECWAM 
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
    SUBROUTINE POSTINTRHS (FL3, FL, IJS, IJL, IG, &
     &                    THWOLD, USOLD, &
     &                    TAUW, Z0OLD, &
     &                    ROAIRO, ZIDLOLD, &
     &                    U10NEW, THWNEW, USNEW, &
     &                    Z0NEW, ROAIRN, ZIDLNEW, &
     &                    SL, FCONST, FMEANWS, MIJ)

       USE DATAPOOL, ONLY : MNP, FR, WETAIL, FRTAIL, WP1TAIL, ISHALLO, FRINTF, COFRM4, CG, WK, &
     &                      DFIM, DFIMOFR, DFFR, DFFR2, WK, RKIND, EMEAN, FMEAN, TH, RKIND, DELU, &
     &                      JUMAX, DT4S, FRM5, IPHYS, LOUTWAM, CD, UFRIC, ALPHA_CH, Z0, ITEST, LCFLX, &
     &                      DELTH => DDIR, TESTNODE, &
     &                      G => G9, &
     &                      ZPI => PI2, &
     &                      EPSMIN => SMALL, &
     &                      NANG => MDC, &
     &                      NFRE => MSC, &
     &                      INDEP => DEP, &
     &                      ROWATER => RHOW, &
     &                      RHOAIR => RHOA

      IMPLICIT NONE

! ----------------------------------------------------------------------
!     ALLOCATED ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      REAL(rkind) :: FL(IJS:IJL,NANG,NFRE),FL3(IJS:IJL,NANG,NFRE),SL(IJS:IJL,NANG,NFRE)
      REAL(rkind),DIMENSION(IJS:IJL,NFRE) :: FCONST
      REAL(rkind),DIMENSION(IJS:IJL) :: THWOLD,USOLD,Z0OLD,TAUW, &
     &                           ROAIRO,ZIDLOLD,FMEANWS
      REAL(rkind),DIMENSION(IJS:IJL) :: U10NEW,THWNEW,USNEW,Z0NEW, &
     &                           ROAIRN,ZIDLNEW

! ----------------------------------------------------------------------

      INTEGER :: IJ,IJS,IJL,K,L,M,IG,ILEV,IDELT,IU06
      INTEGER :: JU(IJS:IJL)
      INTEGER :: MIJ(IJS:IJL)
      REAL(rkind) :: GTEMP1, GTEMP2, FLHAB, XJ, DELT, DELT5, XIMP, AKM1
      REAL(rkind) :: AK2VGM1, XN, PHIDIAG, TAU
      REAL(rkind), DIMENSION(NFRE) :: DELFL
      REAL(rkind), DIMENSION(IJS:IJL) :: EMEANWS, USFM, GADIAG
      REAL(rkind), DIMENSION(IJS:IJL) :: F1MEAN, AKMEAN, XKMEAN
      REAL(rkind), DIMENSION(IJS:IJL) :: PHIEPS, TAUOC, PHIAW
      REAL(rkind), DIMENSION(IJS:IJL) :: TAUWLF,TAUWD,PHIAWDIAG,PHIAWUNR,PHIOC,PHIWA
      REAL(rkind), DIMENSION(IJS:IJL,NANG) :: SPRD
      REAL(rkind), DIMENSION(IJS:IJL,NFRE) :: TEMP, TEMP2
      REAL(rkind), DIMENSION(IJS:IJL,NANG,NFRE) :: CIWAB
      REAL(rkind), DIMENSION(IJS:IJL,NANG,NFRE) :: XLLWS

      IDELT = INT(DT4S)
      ILEV=1

      DO IJ=IJS,IJL
        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'INIT OF POST SOURCE TERMS '
        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(5F20.10)') SUM(FL3) 
      ENDDO

      IF(ISHALLO.EQ.1) THEN
        DO M=1,NFRE
          DO IJ=IJS,IJL
            TEMP2(IJ,M) = FRM5(M)
          ENDDO
        ENDDO
      ELSE
        DO M=1,NFRE
          DO IJ=IJS,IJL
!AR: WAM TABLE REPLACES BY WWM WK            AKM1 = 1./TFAK(INDEP(IJ),M)
            AKM1 = 1./WK(IJ,M)
!AR: WAM TABLE REPLACES BY WWM CG            AK2VGM1 = AKM1**2/TCGOND(INDEP(IJ),M)
            AK2VGM1 = AKM1**2/CG(IJ,M)
            TEMP2(IJ,M) = AKM1*AK2VGM1
!            WRITE(111113,'(4F20.10)') AKM1, AK2VGM1, TEMP2(IJ,M) 
          ENDDO
        ENDDO
      ENDIF

      IF(IPHYS.EQ.0) THEN
        CALL SINPUT (FL3, FL, IJS, IJL, THWNEW, USNEW, Z0NEW, &
     &             ROAIRN, ZIDLNEW, SL, XLLWS)
      ELSE
        CALL SINPUT_ARD (FL3, FL, IJS, IJL, THWNEW, USNEW, Z0NEW, &
     &             ROAIRN, ZIDLNEW, SL, XLLWS)
      ENDIF

!*    2.5 REPLACE DIAGNOSTIC PART OF SPECTRA BY A F**(-5) TAIL.
!         -----------------------------------------------------


!*    2.5.1 COMPUTE MEAN PARAMETERS.
!           ------------------------

      CALL FKMEAN(FL3, IJS, IJL, EMEAN(IJS), FMEAN(IJS), &
     &            F1MEAN, AKMEAN, XKMEAN)

!     MEAN FREQUENCY CHARACTERISTIC FOR WIND SEA
      CALL FEMEANWS(FL3,IJS,IJL,EMEANWS,FMEANWS,XLLWS)

!*    2.5.3 COMPUTE TAIL ENERGY RATIOS.
!           ---------------------------

!     COMPUTE LAST FREQUENCY INDEX OF PROGNOSTIC PART OF SPECTRUM.
      CALL FRCUTINDEX(IJS, IJL, FMEAN(IJS), FMEANWS, MIJ)

      DO IJ=IJS,IJL
        GADIAG(IJ) = 1./TEMP2(IJ,MIJ(IJ))
        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER MEAN PARAMETER'
        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,5F20.10)') MIJ(IJ), AKMEAN, FMEANWS, TEMP2(IJ,MIJ(IJ)), GADIAG(IJ)
      ENDDO


!*    2.5.4 MERGE TAIL INTO SPECTRA.
!           ------------------------

      DO IJ=IJS,IJL
        DO M=1,MIJ(IJ)
          FCONST(IJ,M) = 1.
          TEMP(IJ,M) = 0.
      !    WRITE(111113,'(I10,2F15.10)') M, FCONST(IJ,M), TEMP(IJ,M)
        ENDDO
        DO M=MIJ(IJ)+1,NFRE
          FCONST(IJ,M) = 0.
          TEMP(IJ,M) = TEMP2(IJ,M)*GADIAG(IJ)
      !WRITE(111113,'(I10,3F15.10)')M,FCONST(IJ,M),TEMP(IJ,M),GADIAG(IJ)
        ENDDO
      ENDDO


      DO K=1,NANG
        DO IJ=IJS,IJL
          GADIAG(IJ) = FL3(IJ,K,MIJ(IJ))
      !    WRITE(111113,'(I10,10F20.10)') K, FL3(IJ,K,MIJ(IJ))
        ENDDO
        DO M=1,NFRE
          DO IJ=IJS,IJL
!AR: ICE            FLLOWEST = FLMINFR(JU(IJ),M)*SPRD(IJ,K)
!AR: ICE            FL3(IJ,K,M) = GADIAG(IJ)*TEMP(IJ,M) &
!AR: ICE     &       + MAX(FL3(IJ,K,M),FLLOWEST)*FCONST(IJ,M)
            FL3(IJ,K,M) = GADIAG(IJ)*TEMP(IJ,M) &
     &       + FL3(IJ,K,M)*FCONST(IJ,M)
!            WRITE(111113,'(2I10,10F20.10)')K,M,FL3(IJ,K,M),&
!     &                   TEMP(IJ,M),GADIAG(IJ),FCONST(IJ,M)
          ENDDO
        ENDDO
      ENDDO

      DO IJ=IJS,IJL
        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'BEFORE WIND INPUT 2'
        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(8F20.10)') SUM(FL3) , SUM(FL), THWNEW(IJ)
        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(8F20.10)') USNEW(IJ), Z0NEW(IJ), ZIDLNEW(IJ)
        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(8F20.10)') SUM(SL), SUM(XLLWS(IJ,:,:))
      ENDDO

!*    2.5.5 REEVALUATE WIND INPUT SOURCE TERM, AND WAVE DISSIPATION.
!           -------------------------------------------------------

      IF(IPHYS.EQ.0) THEN
        CALL SINPUT (FL3, FL, IJS, IJL, THWNEW, USNEW, Z0NEW, &
     &             ROAIRN, ZIDLNEW, SL, XLLWS) 
      ELSE
        CALL SINPUT_ARD (FL3, FL, IJS, IJL, THWNEW, USNEW, Z0NEW, &
     &             ROAIRN, ZIDLNEW, SL, XLLWS)
      ENDIF

      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER SINPUT 2'
      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(SL)

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SINPUT CALLED AT THE END'
        CALL FLUSH (IU06)
      ENDIF
      CALL STRESSO (FL3, IJS, IJL, THWNEW, USNEW, Z0NEW, &
     &              ROAIRN, TAUW, TAUWLF, PHIWA, &
     &              PHIAWDIAG, PHIAWUNR, SL, MIJ, &
     &              LCFLX)

      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER STRESSO 2'
      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,15F15.7)') IJS, IJL, SUM(FL3), &
     &              THWNEW, USNEW, Z0NEW, &
     &              ROAIRN, TAUW, TAUWLF, PHIWA, &
     &              PHIAWDIAG, PHIAWUNR, SUM(SL), &
     &              MIJ(IJS:IJL)


      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: STRESSO CALLED AT THE END'
        CALL FLUSH (IU06)
      ENDIF

      CALL AIRSEA (U10NEW(IJS), TAUW(IJS), USNEW(IJS), Z0NEW(IJS), &
     & IJS, IJL, ILEV)

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: AIRSEA CALLED AT THE END'
        CALL FLUSH (IU06)
      ENDIF

      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER AIRSEA 3'
      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,10F15.7)') IJS, U10NEW(IJS), TAUW(IJS), &
      &                              USNEW(IJS), Z0NEW(IJS), ILEV

      IF(IPHYS.EQ.0) THEN
        CALL SDISSIP (FL3 ,FL, IJS, IJL, IG, SL, F1MEAN, XKMEAN, &
     &                PHIOC, TAUWD, MIJ)
      ELSE
        CALL SDISS_ARDH_VEC (FL3 ,FL, IJS, IJL, SL, F1MEAN, XKMEAN, &
     &                PHIOC, TAUWD, MIJ)
      ENDIF
      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'AFTER DISSIP' 
      IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(I10,10F15.7)') IJS, SUM(FL), SUM(FL3), SUM(SL) 
      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. IMPLSCH: SDISSIP CALLED AT THE END'
        CALL FLUSH (IU06)
      ENDIF

!*    2.5.6 DETERMINE FLUXES FROM AIR TO WAVE AND FROM WAVE TO OCEAN.
!           -------------------------------------------------------

      IF(LCFLX) THEN
        DO IJ=IJS,IJL
          TAU       = ROAIRN(IJ)*USNEW(IJ)**2
          XN        = ROAIRN(IJ)*USNEW(IJ)**3

          PHIDIAG    = PHIAWDIAG(IJ)+PHIAWUNR(IJ)
          PHIEPS(IJ) = (PHIOC(IJ)-PHIDIAG)/XN 
          PHIAW(IJ)  = (PHIWA(IJ)+PHIAWUNR(IJ))/XN
          TAUOC(IJ)  = (TAU-TAUWLF(IJ)-TAUWD(IJ))/TAU
        ENDDO
      ENDIF

! ----------------------------------------------------------------------

!*    2.6 SAVE WINDS INTO INTERMEDIATE STORAGE.
!         -------------------------------------

      DO IJ=IJS,IJL
        USOLD(IJ) = USNEW(IJ)
        Z0OLD(IJ) = Z0NEW(IJ)
        ROAIRO(IJ) = ROAIRN(IJ)
        ZIDLOLD(IJ) = ZIDLNEW(IJ)
      ENDDO

      DO IJ=IJS,IJL
        UFRIC(IJ) = USNEW(IJ)
        Z0(IJ)    = Z0NEW(IJ)
        CD(IJ)    = (USNEW(IJ)/U10NEW(IJ))**2
        ALPHA_CH(IJ) = G*Z0NEW(IJ)/USNEW(IJ)**2
      ENDDO

! ----------------------------------------------------------------------

      !IF (LHOOK) CALL DR_HOOK('IMPLSCH',1,ZHOOK_HANDLE)

      DO IJ=IJS,IJL
        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,*) 'END OF POST SOURCE TERMS '
        IF (LOUTWAM .AND. IJS == TESTNODE) WRITE(111113,'(5F20.10)') SUM(FL3), UFRIC(IJ), Z0(IJ), CD(IJ)
      ENDDO

      RETURN
      END SUBROUTINE POSTINTRHS 
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
