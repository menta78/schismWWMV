      SUBROUTINE sed_settleveloc(ised)
!--------------------------------------------------------------------!
! This subroutine computes settling velocity from Soulsby (1997)     !
!                                                                    !
! Author: florian ganthy (fganthy@lnec.pt ; florian.ganthy@gmail.com)!
! Date:   21/12/2012                                                 !
!                                                                    !
! History:                                                           !
!                                                                    !
!--------------------------------------------------------------------!

      USE elfe_glbl, ONLY : rkind
      USE sed_mod,   ONLY : rhom,g,Srho,Sd50,Wsed

      IMPLICIT NONE

!- Local variables --------------------------------------------------!
      
      INTEGER,INTENT(IN)    :: ised

      REAL(rkind),PARAMETER :: nu=1.36d-6 ! Cinematic viscosity
      REAL(rkind)           :: ratio,dstar

!- Start Statement --------------------------------------------------!

! - Preliminary parameters
      ratio    = Srho(ised)/rhom
      dstar    = Sd50(ised) * (g*(ratio-1.d0)/nu**2.d0)**(1.d0/3d0)

! - Settling velocity (in m.s-1)
      Wsed(ised) = (nu/Sd50(ised)) * ((10.36d0**2.d0 + 1.049d0*      &
      &                                dstar**3.d0)**0.5d0 - 10.36d0)

!--------------------------------------------------------------------!
      END SUBROUTINE sed_settleveloc    
