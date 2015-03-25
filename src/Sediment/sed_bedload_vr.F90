      SUBROUTINE sed_bedload_vr(ised,inea,dave)
!--------------------------------------------------------------------!
! This routine computes bedload according to van Rijn                !
!                                                                    !
! Author: Knut Kramer                                                !
! Date: 27/11/2012                                                   !
!                                                                    !
! History: 2012/12 - F.Ganthy : form homogenisation of sediments     !
!          routines                                                  !
!                                                                    !
!--------------------------------------------------------------------!

      USE sed_mod
      USE elfe_glbl, ONLY: rkind,errmsg,dpe,nea,dt,nm,eta2
      USE elfe_msgp, ONLY: myrank,parallel_abort

      IMPLICIT NONE

!- Local variables --------------------------------------------------!

      INTEGER,INTENT(IN)  :: ised,inea
      ! depth averaged elementwise hvel
      REAL(rkind),INTENT(IN) :: dave(nea)

      ! arbitrary coefficients
      REAL(rkind) :: cff, cff1, cff2
      REAL(rkind) :: a_slopex, a_slopey
      ! Total water depth
      REAL(rkind) :: htot
      ! Param for long bed slope effects
      REAL(rkind), PARAMETER :: alpha_bs = 1.d0
      ! Param for orth bed slope effects
      REAL(rkind), PARAMETER :: alpha_bn = 1.5d0

!- Start Statement --------------------------------------------------!

      ucr  = 0.d0
      cff  = 0.d0
      cff1 = 0.d0
      cff2 = 0.d0

!---------------------------------------------------------------------
! - Critical velocity for currents based on Shields (initiation of 
! motion), (van Rijn 2007a)
! Ucr = 0.19*d50**0.1*log10(4*h/12*d90) if 5e-5 < d50 < 5e-4
! Ucr = 8.5*d50**0.6*log10(4*h/12*d90) if 5e-4<= d50 <2e-3
! We do not have d90 in current implementation
!---------------------------------------------------------------------

      htot = dpe(inea)+(eta2(nm(inea,1))+   &
      &                 eta2(nm(inea,2))+   &
      &                 eta2(nm(inea,3)))/3.0d0

      IF (Sd50(ised)>5.0d-5.AND.Sd50(ised)<5.0d-4) THEN

        ucr = 0.19d0*Sd50(ised)**0.1d0*LOG10(4.0d0*htot/Sd50(ised))

      ELSEIF (Sd50(ised)>=5.0d-4.AND.Sd50(ised)<2.0d-3)THEN

        ucr = 8.5d0*Sd50(ised)**0.6d0*LOG10(4.0d0*htot/Sd50(ised))

      ELSE

        WRITE(errmsg,*)'Sediment diameter out of range:',ised,       &
        &              Sd50(ised)
        CALL parallel_abort(errmsg)

      ENDIF

!---------------------------------------------------------------------
! - Mobility parameter (van Rijn 2007a)
! All terms related to waves stripped out...
!---------------------------------------------------------------------

      me = (dave(inea)-ucr)/SQRT(smgd)

!---------------------------------------------------------------------
! - Simplified bed load transport formula for steady flow 
! (van Rijn, 2007)
! rho_s (van Rijn) missing at this point
! instead of [kg/m/s] (van Rijn) here [m2/s]
!
!jl. Looks to me like the units are [m2/s] 
!    In ROMS the bedld is integrated in time(s) and space(1-d, size of
!    RHO-grid spacing to end up with [kg]. Here it looks like rho is 
!    not used yet so the JCG can be used to calculate the change 
!    layer thickness due bedload flux.
!
!    rho is eventually considered at 762
!
! No wave terms implemented. 
!---------------------------------------------------------------------

     bedld = MAX(0.015d0*dave(inea)*htot * &
     &           (Sd50(ised)/htot)**1.2d0*me**1.5d0, 0.0d0)

!---------------------------------------------------------------------
! - Partition bedld into x and y directions, at the center of each
! element (rho points), and integrate in time.
! FX_r and FY_r have dimensions of [m2]
!---------------------------------------------------------------------

      FX_r(inea) = bedld*angleu*dt
      FY_r(inea) = bedld*anglev*dt


!---------------------------------------------------------------------
! why check for isnan == -1?
! what does it mean -1?
!
!jl.  
! This is a non-standard implementation of nan and is compiler dependent.
! Standard IEEE_isnan is only included in 2003+.
! gfortran isnan returns logical
! digital fortran returns logical
! ifort returns logical...
! Compaq visual fortran .true. == -1, so perhaps that was the reasoning. 
! All ==-1 were changed to ==.true. 
!---------------------------------------------------------------------

      IF(isnan(FX_r(inea)).EQV..TRUE.) THEN
        WRITE(errmsg,*)'FX_r0 is NaN',myrank,inea,FX_r(inea),bedld,  &
        &              angleu,dt
        CALL parallel_abort(errmsg)
      ENDIF

      IF(isnan(FY_r(inea)).EQV..TRUE.) THEN
        WRITE(errmsg,*)'FY_r0 is NaN',myrank,inea,FY_r(inea),bedld,  &
        &              anglev,dt
        CALL parallel_abort(errmsg)
      endif

!---------------------------------------------------------------------
! - Bed_slope effects
! longitudinal bed slope
! limit slope to 0.9*(sed_angle)
!---------------------------------------------------------------------

      cff  = (dzdx*angleu+dzdy*anglev)
      cff1 = MIN(ABS(cff),0.9d0*sed_angle)*SIGN(1.0d0,cff)
      cff2 = DATAN(cff1)
      a_slopex = 1.0d0+alpha_bs*                                     &
      &          ((sed_angle/(COS(cff2)*(sed_angle-cff1)))-1.0d0)

!---------------------------------------------------------------------
! - Add contribution of longitudal bed slope to bedload transport
!---------------------------------------------------------------------

      FX_r(inea) = FX_r(inea)*a_slopex
      FY_r(inea) = FY_r(inea)*a_slopex

!---------------------------------------------------------------------
! - Transverse bed slope
!---------------------------------------------------------------------

      cff = (-(dzdx*anglev)+dzdy*angleu)
      cff1 = ABS(bustr(inea))+ABS(bvstr(inea))

      ! - Test used to prevent Inf & NaNs with very small bustr/bvstr 
      !   May still produce unrealistic values 
      IF(cff1<1d-10) THEN
        a_slopey = 0.d0
      ELSE
        cff2 = SQRT(tau_ce(ised)/cff1)
        a_slopey=alpha_bn*cff2*cff
      ENDIF

!---------------------------------------------------------------------
! - Add contribution of transverse to bed load 
!---------------------------------------------------------------------

      FX_r(inea) = FX_r(inea)-(FY_r(inea)*a_slopey)
      FY_r(inea) = FY_r(inea)+(FX_r(inea)*a_slopey)

!---------------------------------------------------------------------
! - Consistency check
!---------------------------------------------------------------------

      IF (isnan(FX_r(inea)).EQV..TRUE.) THEN
        WRITE(errmsg,'(A,I5,A,E15.8,A,E15.8,A,E15.8,A,E15.8,A,E15.8)')&
        &     'FX_r NaN nea:',inea,' bustr:',bustr(inea),             &
        &' bvstr:',bvstr(inea),' cff:',cff,' cff1:',cff1,' cff2:',cff2
        CALL parallel_abort(errmsg)
      ENDIF

      IF (isnan(FY_r(inea)).EQV..TRUE.) THEN
        WRITE(errmsg,*)'FY_r1 is NaN',myrank,inea,FY_r(inea),        &
        &              a_slopey,a_slopex
        CALL parallel_abort(errmsg)
      ENDIF

!---------------------------------------------------------------------
      END SUBROUTINE sed_bedload_vr
