      ! author: Lorenzo Mentaschi
      !  Based on UOST (Mentaschi et al. 2015, 2018, 2020)
      !  
      !  From the ice concentration an isotropic transparency coefficient is estimated.
      !  To simplify the problem, for now it is assumed that beta==alpha
      !  which means that all the energy is dissipated in the current cell.
      !  This gets rid of the shadow problems.
      !  This is inaccurate for the current cell, but should dissipate all the energy needed.
      !  Furthermore, local wave growth (which reduces the effect of unresolved obstacles) 
      !  is neglected
      !
      !  In the future, if this approximation will not be enough, 
      !  some improvement could be introduced:
      !  - a beta different from alpha could be estimated by assuming 
      !  a uniform distribution of the ice in the cell
      !  and the shadow could be estimated for the neighboring cells.
      !  - local wave growth could be taken into account (see the psi function in UOST)
      !  - if there is information on the size of the ice flows, the transparency coeff. 
      !  could be made frequency-dependent.
      SUBROUTINE ICEDISSIP_SRCTRM(IP, SPEC, S, D)
         USE DATAPOOL

         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP ! this is the node id in the local partition
         REAL(rkind), INTENT(IN) :: SPEC(NUMSIG, NUMDIR)
         REAL(rkind), INTENT(INOUT) :: S(NUMSIG, NUMDIR), D(NUMSIG, NUMDIR)
         REAL(rkind) :: ICEC, FREESURF, BETA, CELLAREA, CELLSIZE, CGI, GAM
         REAL(rkind) :: GAMMAUP = 20
         INTEGER  :: IK
   
         S = 0
         D = 0
         ICEC = ICECONC(IP)
         ! computing the transparency coefficient
         ! the total transparency coefficient alpha is given by sqrt(1-concentration)
         ! here we assume that beta==alpha
         FREESURF = MAX(MIN(1-ICEC, 1.), 0.)
         BETA = SQRT(FREESURF)

         IF (BETA .EQ. 1) RETURN

         CELLAREA = SI(IP)
         ! cellsize computed as the diameter of the equivalent circle
         CELLSIZE = SQRT(4/PI*CELLAREA) 

         GAM = (1 - BETA)/BETA
         GAM = MIN(GAM, GAMMAUP)
   
         DO IK = 1,NUMSIG
           CGI = CG(IK,IP)
           D(IK, :) = - CGI/CELLSIZE * GAM
         END DO
         S = D*SPEC
         
      END SUBROUTINE ICEDISSIP_SRCTRM

     
