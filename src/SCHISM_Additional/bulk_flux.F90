      MODULE bulk_flux_mod
        LOGICAL L_COOL_SKIN = .TRUE.
        LOGICAL L_COARE_OOST = .FALSE.
        LOGICAL L_COARE_TAYLOR_YELLAND = .FALSE.
        LOGICAL L_LONGWAVE=.FALSE.
        LOGICAL L_LONGWAVE_OUT=.TRUE.
!
!svn $Id: bulk_flux.F 751 2015-01-07 22:56:36Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2015 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the bulk parameterization of surface wind     !
!  stress and surface net heat fluxes.                                 !
!                                                                      !
!  References:                                                         !
!                                                                      !
!    Fairall, C.W., E.F. Bradley, D.P. Rogers, J.B. Edson and G.S.     !
!      Young, 1996:  Bulk parameterization of air-sea fluxes for       !
!      tropical ocean-global atmosphere Coupled-Ocean Atmosphere       !
!      Response Experiment, JGR, 101, 3747-3764.                       !
!                                                                      !
!    Fairall, C.W., E.F. Bradley, J.S. Godfrey, G.A. Wick, J.B.        !
!      Edson, and G.S. Young, 1996:  Cool-skin and warm-layer          !
!      effects on sea surface temperature, JGR, 101, 1295-1308.        !
!                                                                      !
!    Liu, W.T., K.B. Katsaros, and J.A. Businger, 1979:  Bulk          !
!        parameterization of the air-sea exchange of heat and          !
!        water vapor including the molecular constraints at            !
!        the interface, J. Atmos. Sci, 36, 1722-1735.                  !
!                                                                      !
!  Adapted from COARE code written originally by David Rutgers and     !
!  Frank Bradley.                                                      !
!                                                                      !
!  EMINUSP option for equivalent salt fluxes added by Paul Goodman     !
!  (10/2004).                                                          !
!                                                                      !
!  Modified by Kate Hedstrom for COARE version 3.0 (03/2005).          !
!  Modified by Jim Edson to correct specific hunidities.               !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!     Fairall et al., 2003: J. Climate, 16, 571-591.                   !
!                                                                      !
!     Taylor, P. K., and M. A. Yelland, 2001: The dependence of sea    !
!     surface roughness on the height and steepness of the waves.      !
!     J. Phys. Oceanogr., 31, 572-590.                                 !
!                                                                      !
!     Oost, W. A., G. J. Komen, C. M. J. Jacobs, and C. van Oort, 2002:!
!     New evidence for a relation between wind stress and wave age     !
!     from measurements during ASGAMAGE. Bound.-Layer Meteor., 103,    !
!     409-438.                                                         !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: bulk_flux, bulk_psiu, bulk_psit
!
      CONTAINS
      SUBROUTINE bulk_flux      (prho, tr_nd,                           &
     &                           Hair_spec, Pair, Tair, Uwind, Vwind,   &
     &                           cloud,                                 &
     &                           rain, lhflx, lrflx, shflx,             &
     &                           srflx, stflx,                          &
#ifdef PREC_EVAP
     &                           EminusP, evap,                         &
#endif
     &                           sustr, svstr)
!
!  Imported variable declarations.
!
      real(rkind), intent(in) :: alpha(npa)
      real(rkind), intent(in) :: beta(npa)
      real(rkind), intent(in) :: prho(nvrt,npa)
      real(rkind), intent(in) :: tr_nd(ntracers,nvrt,npa)
      real(rkind), intent(in) :: Hair_spec(npa)
      real(rkind), intent(in) :: Pair(npa)
      real(rkind), intent(in) :: Tair(npa)
      real(rkind), intent(in) :: Uwind(npa)
      real(rkind), intent(in) :: Vwind(npa)
      real(rkind), intent(in) :: cloud(npa)
      real(rkind), intent(in) :: rain(npa)

      real(rkind), intent(inout) :: lhflx(npa)
      real(rkind), intent(inout) :: lrflx(npa)
      real(rkind), intent(inout) :: shflx(npa)
      real(rkind), intent(inout) :: srflx(npa)
      real(rkind), intent(inout) :: stflx(npa,NT(ng))

      real(rkind), intent(out) :: EminusP(npa)
      real(rkind), intent(out) :: evap(npa)
      real(rkind), intent(out) :: sustr(npa)
      real(rkind), intent(out) :: svstr(npa)
!
!  Local variable declarations.
!
      integer :: Iter, i, j, k
      integer, parameter :: IterMax = 3

      real(rkind), parameter :: eps = 1.0E-20_rkind
      real(rkind), parameter :: r3 = 1.0_rkind/3.0_rkind

      real(rkind), parameter :: blk_Cpa = 1004.67_rkind      ! (J/kg/K), Businger 1982
      real(rkind), parameter :: blk_Cpw = 4000.0_rkind       ! (J/kg/K)
      real(rkind), parameter :: blk_Rgas = 287.1_rkind       ! (J/kg/K)
      real(rkind), parameter :: blk_Zabl = 600.0_rkind       ! (m)
      real(rkind), parameter :: blk_beta = 1.2_rkind         ! non-dimensional
      real(rkind), parameter :: blk_dter = 0.3_rkind         ! (K)
      real(rkind), parameter :: blk_tcw = 0.6_rkind          ! (W/m/K)
      real(rkind), parameter :: blk_visw = 0.000001_rkind    ! (m2/s)
      real(rkind), parameter :: Cp = 3985.0_rkind              ! Joules/kg/degC
      real(rkind), parameter :: Csolar = 1353.0_rkind          ! 1360-1380 W/m2
      real(rkind), parameter :: Eradius = 6371315.0_rkind      ! m
      real(rkind), parameter :: StefBo = 5.67E-8_rkind         ! Watts/m2/K4
      real(rkind), parameter :: emmiss = 0.97_rkind            ! non_dimensional
      real(rkind), parameter :: rhow = 1000.0_rkind            ! kg/m3
      real(rkind), parameter :: g = 9.81_rkind                 ! m/s2
      real(rkind), parameter :: gorho0                      ! m4/s2/kg
      real(rkind), parameter :: vonKar = 0.41_rkind            ! non-dimensional
      

      
! alpha and beta are actualla variable. But here we use fixed value
      real(rkind), parameter :: alpha = 2.1014611551470d-04
      real(rkind), parameter :: beta = 7.2575037309946d-04

      real(rkind) :: Bf, Cd, Hl, Hlw, Hscale, Hs, Hsr, IER
      real(rkind) :: PairM,  RH, Taur
      real(rkind) :: Wspeed, ZQoL, ZToL

      real(rkind) :: cff, cff1, cff2, diffh, diffw, oL, upvel
      real(rkind) :: twopi_inv, wet_bulb
      real(rkind) :: e_sat, vap_p
      real(rkind) :: Clam, Fc, Hcool, Hsb, Hlb, Qbouy, Qcool, lambd

      real(rkind), dimension(npa) :: CC
      real(rkind), dimension(npa) :: Cd10
      real(rkind), dimension(npa) :: Ch10
      real(rkind), dimension(npa) :: Ct10
      real(rkind), dimension(npa) :: charn
      real(rkind), dimension(npa) :: Ct
      real(rkind), dimension(npa) :: Cwave
      real(rkind), dimension(npa) :: Cwet
      real(rkind), dimension(npa) :: delQ
      real(rkind), dimension(npa) :: delQc
      real(rkind), dimension(npa) :: delT
      real(rkind), dimension(npa) :: delTc
      real(rkind), dimension(npa) :: delW
      real(rkind), dimension(npa) :: L
      real(rkind), dimension(npa) :: L10
      real(rkind), dimension(npa) :: Q
      real(rkind), dimension(npa) :: Qair
      real(rkind), dimension(npa) :: Qpsi
      real(rkind), dimension(npa) :: Qsea
      real(rkind), dimension(npa) :: Qstar
      real(rkind), dimension(npa) :: rhoAir
      real(rkind), dimension(npa) :: rhoSea
      real(rkind), dimension(npa) :: Ri
      real(rkind), dimension(npa) :: Ribcu
      real(rkind), dimension(npa) :: Rr
      real(rkind), dimension(npa) :: Scff
      real(rkind), dimension(npa) :: TairC
      real(rkind), dimension(npa) :: TairK
      real(rkind), dimension(npa) :: Tcff
      real(rkind), dimension(npa) :: Tpsi
      real(rkind), dimension(npa) :: TseaC
      real(rkind), dimension(npa) :: TseaK
      real(rkind), dimension(npa) :: Tstar
      real(rkind), dimension(npa) :: u10
      real(rkind), dimension(npa) :: VisAir
      real(rkind), dimension(npa) :: WaveLength
      real(rkind), dimension(npa) :: Wgus
      real(rkind), dimension(npa) :: Wmag
      real(rkind), dimension(npa) :: Wpsi
      real(rkind), dimension(npa) :: Wstar
      real(rkind), dimension(npa) :: Zetu
      real(rkind), dimension(npa) :: Zo10
      real(rkind), dimension(npa) :: ZoT10
      real(rkind), dimension(npa) :: ZoL
      real(rkind), dimension(npa) :: ZoQ
      real(rkind), dimension(npa) :: ZoT
      real(rkind), dimension(npa) :: ZoW
      real(rkind), dimension(npa) :: ZWoL

      real(rkind), dimension(npa) :: Hlv
      real(rkind), dimension(npa) :: LHeat
      real(rkind), dimension(npa) :: LRad
      real(rkind), dimension(npa) :: SHeat
      real(rkind), dimension(npa) :: SRad
      real(rkind), dimension(npa) :: Taux
      real(rkind), dimension(npa) :: Tauy
!
!=======================================================================
!  Atmosphere-Ocean bulk fluxes parameterization.
!=======================================================================
!
      Hscale=rho0*Cp
      twopi_inv=0.5_rkind/pi
      DO i=1,npa
          IF (idry(i) == 1) CYCLE
!
!  Input bulk parameterization fields.
!
          Wmag(i)=SQRT(Uwind(i)*Uwind(i)+Vwind(i)*Vwind(i))
          PairM=Pair(i) / 100.0_rkind
          TairC(i)=Tair(i)
          TairK(i)=TairC(i)+273.16_rkind
          TseaC(i)=tr_nd(1,nvrt,i)
          TseaK(i)=TseaC(i)+273.16_rkind
          rhoSea(i)=prho(nvrt,i)
          SpecHum=Hair_spec(i)
          RH  = -1000000    ! We need to put conversion of specific humidity to
          ! relative humidity but that is needed only for LONGWAVE option
          ! which is questionable
          SRad(i)=srflx(i)*Hscale
          Tcff(i)=alpha
          Scff(i)=beta
!
!  Initialize.
!
          delTc(i)=0.0_rkind
          delQc(i)=0.0_rkind
          LHeat(i)=lhflx(i)*Hscale
          SHeat(i)=shflx(i)*Hscale
          Taur=0.0_rkind
          Taux(i)=0.0_rkind
          Tauy(i)=0.0_rkind
!
!-----------------------------------------------------------------------
!  Compute net longwave radiation (W/m2), LRad.
!-----------------------------------------------------------------------

       IF (L_LONGWAVE) THEN
!
!  Use Berliand (1952) formula to calculate net longwave radiation.
!  The equation for saturation vapor pressure is from Gill (Atmosphere-
!  Ocean Dynamics, pp 606). Here the coefficient in the cloud term
!  is assumed constant, but it is a function of latitude varying from
!  1.0 at poles to 0.5 at the Equator).
!
          cff=(0.7859_rkind+0.03477_rkind*TairC(i))/                          &
     &        (1.0_rkind+0.00412_rkind*TairC(i))
          e_sat=10.0_rkind**cff   ! saturation vapor pressure (hPa or mbar)
          vap_p=e_sat*RH       ! water vapor pressure (hPa or mbar)
          cff2=TairK(i)*TairK(i)*TairK(i)
          cff1=cff2*TairK(i)
          LRad(i)=-emmiss*StefBo*                                     &
     &              (cff1*(0.39_rkind-0.05_rkind*SQRT(vap_p))*                &
     &                    (1.0_rkind-0.6823_rkind*cloud(i)*cloud(i))+     &
     &               cff2*4.0_rkind*(TseaK(i)-TairK(i)))

        ELSE IF (L_LONGWAVE_OUT) THEN
!
!  Treat input longwave data as downwelling radiation only and add
!  outgoing IR from model sea surface temperature.
!
          LRad(i)=lrflx(i)*Hscale-                                  &
     &              emmiss*StefBo*TseaK(i)*TseaK(i)*TseaK(i)*TseaK(i)

        ELSE
          LRad(i)=lrflx(i)*Hscale
        END IF
!
!-----------------------------------------------------------------------
!  Compute specific humidities (kg/kg).
!
!    note that Qair(i) is the saturation specific humidity at Tair
!                 Q(i) is the actual specific humidity
!              Qsea(i) is the saturation specific humidity at Tsea
!
!          Saturation vapor pressure in mb is first computed and then
!          converted to specific humidity in kg/kg
!
!          The saturation vapor pressure is computed from Teten formula
!          using the approach of Buck (1981):
!
!          Esat(mb) = (1.0007_rkind+3.46E-6_rkind*PairM(mb))*6.1121_rkind*
!                  EXP(17.502_rkind*TairC(C)/(240.97_rkind+TairC(C)))
!
!          The ambient vapor is found from the definition of the
!          Relative humidity:
!
!          RH = W/Ws*100 ~ E/Esat*100   E = RH/100*Esat if RH is in %
!                                       E = RH*Esat     if RH fractional
!
!          The specific humidity is then found using the relationship:
!
!          Q = 0.622 E/(P + (0.622-1)e)
!
!          Q(kg/kg) = 0.62197_rkind*(E(mb)/(PairM(mb)-0.378_rkind*E(mb)))
!
!-----------------------------------------------------------------------
!
!  Compute air saturation vapor pressure (mb), using Teten formula.
!
          cff=(1.0007_rkind+3.46E-6_rkind*PairM)*6.1121_rkind*                   &
     &        EXP(17.502_rkind*TairC(i)/(240.97_rkind+TairC(i)))
!
!  Compute specific humidity at Saturation, Qair (kg/kg).
!
          Qair(i)=0.62197_rkind*(cff/(PairM-0.378_rkind*cff))
!
!  Compute specific humidity, Q (kg/kg).
!
          Q(i)=SpecHum
!          IF (RH.lt.2.0_rkind) THEN                       !RH fraction
!            cff=cff*RH                                 !Vapor pres (mb)
!            Q(i)=0.62197_rkind*(cff/(PairM-0.378_rkind*cff)) !Spec hum (kg/kg)
!          ELSE          !RH input was actually specific humidity in g/kg
!            Q(i)=RH/1000.0_rkind                          !Spec Hum (kg/kg)
!          END IF
!
!  Compute water saturation vapor pressure (mb), using Teten formula.
!
          cff=(1.0007_rkind+3.46E-6_rkind*PairM)*6.1121_rkind*                   &
     &        EXP(17.502_rkind*TseaC(i)/(240.97_rkind+TseaC(i)))
!
!  Vapor Pressure reduced for salinity (Kraus & Businger, 1994, pp 42).
!
          cff=cff*0.98_rkind
!
!  Compute Qsea (kg/kg) from vapor pressure.
!
          Qsea(i)=0.62197_rkind*(cff/(PairM-0.378_rkind*cff))
!
!-----------------------------------------------------------------------
!  Compute Monin-Obukhov similarity parameters for wind (Wstar),
!  heat (Tstar), and moisture (Qstar), Liu et al. (1979).
!-----------------------------------------------------------------------
!
!  Moist air density (kg/m3).
!
          rhoAir(i)=PairM*100.0_rkind/(blk_Rgas*TairK(i)*                  &
     &                              (1.0_rkind+0.61_rkind*Q(i)))
!
!  Kinematic viscosity of dry air (m2/s), Andreas (1989).
!
          VisAir(i)=1.326E-5_rkind*                                        &
     &              (1.0_rkind+TairC(i)*(6.542E-3_rkind+TairC(i)*             &
     &               (8.301E-6_rkind-4.84E-9_rkind*TairC(i))))
!
!  Compute latent heat of vaporization (J/kg) at sea surface, Hlv.
!
          Hlv(i)=(2.501_rkind-0.00237_rkind*TseaC(i))*1.0E+6_rkind
!
!  Assume that wind is measured relative to sea surface and include
!  gustiness.
!
          Wgus(i)=0.5_rkind
          delW(i)=SQRT(Wmag(i)*Wmag(i)+Wgus(i)*Wgus(i))
          delQ(i)=Qsea(i)-Q(i)
          delT(i)=TseaC(i)-TairC(i)
!
!  Neutral coefficients.
!
          ZoW(i)=0.0001_rkind
          u10(i)=delW(i)*LOG(10.0_rkind/ZoW(i))/LOG(blk_ZW(ng)/ZoW(i))
          Wstar(i)=0.035_rkind*u10(i)
          Zo10(i)=0.011_rkind*Wstar(i)*Wstar(i)/g+                         &
     &            0.11_rkind*VisAir(i)/Wstar(i)
          Cd10(i)=(vonKar/LOG(10.0_rkind/Zo10(i)))**2
          Ch10(i)=0.00115_rkind
          Ct10(i)=Ch10(i)/sqrt(Cd10(i))
          ZoT10(i)=10.0_rkind/EXP(vonKar/Ct10(i))
          Cd=(vonKar/LOG(blk_ZW(ng)/Zo10(i)))**2
!
!  Compute Richardson number.
!
          Ct(i)=vonKar/LOG(blk_ZT(ng)/ZoT10(i))  ! T transfer coefficient
          CC(i)=vonKar*Ct(i)/Cd
          delTc(i)=0.0_rkind
          IF (L_COOL_SKIN) THEN
            delTc(i)=blk_dter
          END IF
          Ribcu(i)=-blk_ZW(ng)/(blk_Zabl*0.004_rkind*blk_beta**3)
          Ri(i)=-g*blk_ZW(ng)*((delT(i)-delTc(i))+                      &
     &                          0.61_rkind*TairK(i)*delQ(i))/              &
     &          (TairK(i)*delW(i)*delW(i))
          IF (Ri(i).lt.0.0_rkind) THEN
            Zetu(i)=CC(i)*Ri(i)/(1.0_rkind+Ri(i)/Ribcu(i))       ! Unstable
          ELSE
            Zetu(i)=CC(i)*Ri(i)/(1.0_rkind+3.0_rkind*Ri(i)/CC(i))   ! Stable
          END IF
          L10(i)=blk_ZW(ng)/Zetu(i)
!
!  First guesses for Monon-Obukhov similarity scales.
!
          Wstar(i)=delW(i)*vonKar/(LOG(blk_ZW(ng)/Zo10(i))-             &
     &                             bulk_psiu(blk_ZW(ng)/L10(i),pi))
          Tstar(i)=-(delT(i)-delTc(i))*vonKar/                          &
     &             (LOG(blk_ZT(ng)/ZoT10(i))-                           &
     &              bulk_psit(blk_ZT(ng)/L10(i),pi))
          Qstar(i)=-(delQ(i)-delQc(i))*vonKar/                          &
     &             (LOG(blk_ZQ(ng)/ZoT10(i))-                           &
     &              bulk_psit(blk_ZQ(ng)/L10(i),pi))
!
!  Modify Charnock for high wind speeds. The 0.125 factor below is for
!  1.0/(18.0-10.0).
!
          IF (delW(i).gt.18.0_rkind) THEN
            charn(i)=0.018_rkind
          ELSE IF ((10.0_rkind.lt.delW(i)).and.(delW(i).le.18.0_rkind)) THEN
            charn(i)=0.011_rkind+                                          &
     &               0.125_rkind*(0.018_rkind-0.011_rkind)*(delW(i)-10._rkind)
          ELSE
            charn(i)=0.011_rkind
          END IF
      END DO
!
!  Iterate until convergence. It usually converges within 3 iterations.
!
      DO Iter=1,IterMax
        DO i=1,npa
            ZoW(i)=charn(i)*Wstar(i)*Wstar(i)/g+                        &
     &             0.11_rkind*VisAir(i)/(Wstar(i)+eps)
            Rr(i)=ZoW(i)*Wstar(i)/VisAir(i)
!
!  Compute Monin-Obukhov stability parameter, Z/L.
!
            ZoQ(i)=MIN(1.15e-4_rkind,5.5e-5_rkind/Rr(i)**0.6_rkind)
            ZoT(i)=ZoQ(i)
            ZoL(i)=vonKar*g*blk_ZW(ng)*                                 &
     &             (Tstar(i)*(1.0_rkind+0.61_rkind*Q(i))+                     &
     &                        0.61_rkind*TairK(i)*Qstar(i))/               &
     &             (TairK(i)*Wstar(i)*Wstar(i)*                         &
     &              (1.0_rkind+0.61_rkind*Q(i))+eps)
            L(i)=blk_ZW(ng)/(ZoL(i)+eps)
!
!  Evaluate stability functions at Z/L.
!
            Wpsi(i)=bulk_psiu(ZoL(i),pi)
            Tpsi(i)=bulk_psit(blk_ZT(ng)/L(i),pi)
            Qpsi(i)=bulk_psit(blk_ZQ(ng)/L(i),pi)
            IF (L_COOL_SKIN) THEN
            Cwet(i)=0.622_rkind*Hlv(i)*Qsea(i)/                          &
     &              (blk_Rgas*TseaK(i)*TseaK(i))
            delQc(i)=Cwet(i)*delTc(i)
            END IF
!
!  Compute wind scaling parameters, Wstar.
!
            Wstar(i)=MAX(eps,delW(i)*vonKar/                            &
     &               (LOG(blk_ZW(ng)/ZoW(i))-Wpsi(i)))
            Tstar(i)=-(delT(i)-delTc(i))*vonKar/                        &
     &               (LOG(blk_ZT(ng)/ZoT(i))-Tpsi(i))
            Qstar(i)=-(delQ(i)-delQc(i))*vonKar/                        &
     &               (LOG(blk_ZQ(ng)/ZoQ(i))-Qpsi(i))
!
!  Compute gustiness in wind speed.
!
            Bf=-g/TairK(i)*                                             &
     &         Wstar(i)*(Tstar(i)+0.61_rkind*TairK(i)*Qstar(i))
            IF (Bf.gt.0.0_rkind) THEN
              Wgus(i)=blk_beta*(Bf*blk_Zabl)**r3
            ELSE
              Wgus(i)=0.2_rkind
            END IF
            delW(i)=SQRT(Wmag(i)*Wmag(i)+Wgus(i)*Wgus(i))
          IF (L_COOL_SKIN) THEN
!
!-----------------------------------------------------------------------
!  Cool Skin correction.
!-----------------------------------------------------------------------
!
!  Cool skin correction constants. Clam: part of Saunders constant
!  lambda; Cwet: slope of saturation vapor.
!
            Clam=16.0_rkind*g*blk_Cpw*(rhoSea(i)*blk_visw)**3.0_rkind/        &
     &           (blk_tcw*blk_tcw*rhoAir(i)*rhoAir(i))
!
!  Set initial guesses for cool-skin layer thickness (Hcool).
!
            Hcool=0.001_rkind
!
!  Backgound sensible and latent heat.
!
            Hsb=-rhoAir(i)*blk_Cpa*Wstar(i)*Tstar(i)
            Hlb=-rhoAir(i)*Hlv(i)*Wstar(i)*Qstar(i)
!
!  Mean absoption in cool-skin layer.
!
            Fc=0.065_rkind+11.0_rkind*Hcool-                                  &
     &         (1.0_rkind-EXP(-Hcool*1250.0_rkind))*6.6E-5_rkind/Hcool
!
!  Total cooling at the interface.
!
            Qcool=LRad(i)+Hsb+Hlb-SRad(i)*Fc
            Qbouy=Tcff(i)*Qcool+Scff(i)*Hlb*blk_Cpw/Hlv(i)
!
!  Compute temperature and moisture change.
!
            IF ((Qcool.gt.0.0_rkind).and.(Qbouy.gt.0.0_rkind)) THEN
              lambd=6.0_rkind/                                             &
     &              (1.0_rkind+                                            &
     &               (Clam*Qbouy/(Wstar(i)+eps)**4.0_rkind)**0.75_rkind)**r3
              Hcool=lambd*blk_visw/(SQRT(rhoAir(i)/rhoSea(i))*          &
     &                              Wstar(i)+eps)
              delTc(i)=Qcool*Hcool/blk_tcw
            ELSE
              delTc(i)=0.0_rkind
            END IF
            delQc(i)=Cwet(i)*delTc(i)
          END IF
        END DO
!
!-----------------------------------------------------------------------
!  Compute Atmosphere/Ocean fluxes.
!-----------------------------------------------------------------------
!
        DO i=1,npa
!
!  Compute transfer coefficients for momentum (Cd).
!
          Wspeed=SQRT(Wmag(i)*Wmag(i)+Wgus(i)*Wgus(i))
          Cd=Wstar(i)*Wstar(i)/(Wspeed*Wspeed+eps)
!
!  Compute turbulent sensible heat flux (W/m2), Hs.
!
          Hs=-blk_Cpa*rhoAir(i)*Wstar(i)*Tstar(i)
!
!  Compute sensible heat flux (W/m2) due to rainfall (kg/m2/s), Hsr.
!
          diffw=2.11E-5_rkind*(TairK(i)/273.16_rkind)**1.94_rkind
          diffh=0.02411_rkind*(1.0_rkind+TairC(i)*                            &
     &                      (3.309E-3_rkind-1.44E-6_rkind*TairC(i)))/         &
     &          (rhoAir(i)*blk_Cpa)
          cff=Qair(i)*Hlv(i)/(blk_Rgas*TairK(i)*TairK(i))
          wet_bulb=1.0_rkind/(1.0_rkind+0.622_rkind*(cff*Hlv(i)*diffw)/        &
     &                                     (blk_Cpa*diffh))
          Hsr=rain(i)*wet_bulb*blk_Cpw*                               &
     &        ((TseaC(i)-TairC(i))+(Qsea(i)-Q(i))*Hlv(i)/blk_Cpa)
          SHeat(i)=(Hs+Hsr)
!
!  Compute turbulent latent heat flux (W/m2), Hl.
!
          Hl=-Hlv(i)*rhoAir(i)*Wstar(i)*Qstar(i)
!
!  Compute Webb correction (Webb effect) to latent heat flux, Hlw.
!
          upvel=-1.61_rkind*Wstar(i)*Qstar(i)-                             &
     &          (1.0_rkind+1.61_rkind*Q(i))*Wstar(i)*Tstar(i)/TairK(i)
          Hlw=rhoAir(i)*Hlv(i)*upvel*Q(i)
          LHeat(i)=(Hl+Hlw)
!
!  Compute momentum flux (N/m2) due to rainfall (kg/m2/s).
!
          Taur=0.85_rkind*rain(i)*Wmag(i)
!
!  Compute wind stress components (N/m2), Tau.
!
          cff=rhoAir(i)*Cd*Wspeed
          Taux(i)=(cff*Uwind(i)+Taur*SIGN(1.0_rkind,Uwind(i)))
          Tauy(i)=(cff*Vwind(i)+Taur*SIGN(1.0_rkind,Vwind(i)))
        END DO
!
!=======================================================================
!  Compute surface net heat flux and surface wind stress.
!=======================================================================
!
!  Compute kinematic, surface, net heat flux (degC m/s).  Notice that
!  the signs of latent and sensible fluxes are reversed because fluxes
!  calculated from the bulk formulations above are positive out of the
!  ocean.
!
!  For EMINUSP option,  EVAP = LHeat (W/m2) / Hlv (J/kg) = kg/m2/s
!                       PREC = rain = kg/m2/s
!
!  To convert these rates to m/s divide by freshwater density, rhow.
!
!  Note that when the air is undersaturated in water vapor (Q < Qsea)
!  the model will evaporate and LHeat > 0:
!
!                   LHeat positive out of the ocean
!                    evap positive out of the ocean
!
!  Note that if evaporating, the salt flux is positive
!        and if     raining, the salt flux is negative
!
!  Note that fresh water flux is positive out of the ocean and the
!  salt flux (stflx(isalt)) is positive into the ocean. It is converted
!  to (psu m/s) for stflx(isalt) in "set_vbc.F". The E-P value is
!  saved in variable EminusP for I/O purposes.
!
      Hscale=1.0_rkind/(rho0*Cp)
      cff=1.0_rkind/rhow
      DO i=1,npa
          lrflx(i)=LRad(i)*Hscale
          lhflx(i)=-LHeat(i)*Hscale
          shflx(i)=-SHeat(i)*Hscale
          stflx(i,itemp)=srflx(i)+lrflx(i)+                      &
     &                      lhflx(i)+shflx(i)
          evap(i)=LHeat(i)/Hlv(i)
          stflx(i,isalt)=cff*(evap(i)-rain(i))
          EminusP(i)=stflx(i,isalt)
      END DO
!
!  Compute kinematic, surface wind stress (m2/s2).
!
      sustr = Taux
      svstr = Tauy
      END SUBROUTINE bulk_flux

      FUNCTION bulk_psiu (ZoL, pi)
!
!=======================================================================
!                                                                      !
!  This function evaluates the stability function for  wind speed      !
!  by matching Kansas  and free convection forms.  The convective      !
!  form follows Fairall et al. (1996) with profile constants from      !
!  Grachev et al. (2000) BLM.  The  stable  form is from Beljaars      !
!  and Holtslag (1991).                                                !
!                                                                      !
!=======================================================================
!
!  Function result
!
      real(rkind) :: bulk_psiu
!
!  Imported variable declarations.
!
      real(rkind), intent(in) :: ZoL, pi
!
!  Local variable declarations.
!
      real(rkind), parameter :: r3 = 1.0_rkind/3.0_rkind

      real(rkind) :: Fw, cff, psic, psik, x, y
!
!-----------------------------------------------------------------------
!  Compute stability function, PSI.
!-----------------------------------------------------------------------
!
!  Unstable conditions.
!
      IF (ZoL.lt.0.0_rkind) THEN
        x=(1.0_rkind-15.0_rkind*ZoL)**0.25_rkind
        psik=2.0_rkind*LOG(0.5_rkind*(1.0_rkind+x))+                             &
     &       LOG(0.5_rkind*(1.0_rkind+x*x))-                                  &
     &       2.0_rkind*ATAN(x)+0.5_rkind*pi
!
!  For very unstable conditions, use free-convection (Fairall).
!
        cff=SQRT(3.0_rkind)
        y=(1.0_rkind-10.15_rkind*ZoL)**r3
        psic=1.5_rkind*LOG(r3*(1.0_rkind+y+y*y))-                             &
     &       cff*ATAN((1.0_rkind+2.0_rkind*y)/cff)+pi/cff
!
!  Match Kansas and free-convection forms with weighting Fw.
!
        cff=ZoL*ZoL
        Fw=cff/(1.0_rkind+cff)
        bulk_psiu=(1.0_rkind-Fw)*psik+Fw*psic
!
!  Stable conditions.
!
      ELSE
        cff=MIN(50.0_rkind,0.35_rkind*ZoL)
        bulk_psiu=-((1.0_rkind+ZoL)+0.6667_rkind*(ZoL-14.28_rkind)/              &
     &            EXP(cff)+8.525_rkind)
      END IF
      RETURN
      END FUNCTION bulk_psiu

      FUNCTION bulk_psit (ZoL, pi)
!
!=======================================================================
!                                                                      !
!  This function evaluates the  stability function  for moisture and   !
!  heat by matching Kansas and free convection forms. The convective   !
!  form follows Fairall et al. (1996) with  profile  constants  from   !
!  Grachev et al. (2000) BLM.  The stable form is from  Beljaars and   !
!  and Holtslag (1991).                                                !
!
!=======================================================================
!
!  Function result
!
      real(rkind) :: bulk_psit
!
!  Imported variable declarations.
!
      real(rkind), intent(in) :: ZoL, pi
!
!  Local variable declarations.
!
      real(rkind), parameter :: r3 = 1.0_rkind/3.0_rkind

      real(rkind) :: Fw, cff, psic, psik, x, y
!
!-----------------------------------------------------------------------
!  Compute stability function, PSI.
!-----------------------------------------------------------------------
!
!  Unstable conditions.
!
      IF (ZoL.lt.0.0_rkind) THEN
        x=(1.0_rkind-15.0_rkind*ZoL)**0.5_rkind
        psik=2.0_rkind*LOG(0.5_rkind*(1.0_rkind+x))
!
!  For very unstable conditions, use free-convection (Fairall).
!
        cff=SQRT(3.0_rkind)
        y=(1.0_rkind-34.15_rkind*ZoL)**r3
        psic=1.5_rkind*LOG(r3*(1.0_rkind+y+y*y))-                             &
     &       cff*ATAN((1.0_rkind+2.0_rkind*y)/cff)+pi/cff
!
!  Match Kansas and free-convection forms with weighting Fw.
!
        cff=ZoL*ZoL
        Fw=cff/(1.0_rkind+cff)
        bulk_psit=(1.0_rkind-Fw)*psik+Fw*psic
!
!  Stable conditions.
!
      ELSE
        cff=MIN(50.0_rkind,0.35_rkind*ZoL)
        bulk_psit=-((1.0_rkind+2.0_rkind*ZoL)**1.5_rkind+                        &
     &            0.6667_rkind*(ZoL-14.28_rkind)/EXP(cff)+8.525_rkind)
      END IF

      RETURN
      END FUNCTION bulk_psit
      END MODULE
