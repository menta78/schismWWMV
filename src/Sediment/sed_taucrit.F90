      SUBROUTINE sed_taucrit(ised)
!--------------------------------------------------------------------!
! This subroutine computes critical shear stres for erosion from     !
! Soulsby (1997)                                                     !
!                                                                    !
! Author: florian ganthy (fganthy@lnec.pt ; florian.ganthy@gmail.com)!
! Date:   21/12/2012                                                 !
!                                                                    !
! History:                                                           !
!                                                                    !
!--------------------------------------------------------------------!

      USE elfe_glbl, ONLY : rkind
      USE sed_mod,   ONLY : rhom,g,Srho,Sd50,tau_ce

      IMPLICIT NONE

!- Local variables --------------------------------------------------!
      
      INTEGER,INTENT(IN)    :: ised

      REAL(rkind),PARAMETER :: nu=1.36d-6 ! Cinematic viscosity
      REAL(rkind)           :: ratio,dstar,theta_cr

!- Start Statement --------------------------------------------------!

! - Preliminary parameters
      ratio    = Srho(ised)/rhom
      dstar    = Sd50(ised) * (g*(ratio-1.d0)/nu**2.d0)**(1.d0/3d0)
      theta_cr = 0.3d0/(1.d0+1.2d0*dstar) + 0.055d0 *                &
      &         (1.d0-EXP(-0.02d0*dstar))

! - Critical shear stress (N.m-2)
      tau_ce(ised) = theta_cr*g*(Srho(ised)-rhom)*Sd50(ised)

!--------------------------------------------------------------------!
      END SUBROUTINE sed_taucrit
