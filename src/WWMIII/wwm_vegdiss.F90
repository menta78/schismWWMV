!****************************************************************
!
!****************************************************************
      SUBROUTINE VEGDISSIP (IP,SSVEG,DSSVEG,ACLOC,DEPTH,ETOT,SBAR,KBAR) 
        USE DATAPOOL
        IMPLICIT NONE
 
        INTEGER, INTENT(IN)     :: IP

        REAL(rkind),INTENT(IN)  :: KBAR,SBAR,ETOT,DEPTH
        REAL(rkind),INTENT(IN)  :: ACLOC(NUMSIG,NUMDIR)
        REAL(rkind),INTENT(OUT) :: SSVEG(NUMSIG,NUMDIR), DSSVEG(NUMSIG,NUMDIR)

        REAL(rkind) :: BGDISS, KDBAR, SVEGET, VEGDISS
        INTEGER     :: IS, ID

        IF (ETOT .LT. THR) RETURN

        KDBAR  = KBAR * DEPTH
!
        BGDISS = SQRT(TWO/PI)*G9**2*(KBAR/SBAR)**3*SQRT(ETOT)/(THREE*KBAR*(COSH(KDBAR))**3)*NPLANTSPSQM(IP)
!
!        vdrgcoeff_local =  ! drag coefficient / per layer and node 
!        vdiam_local     =  ! diam of veg. / per layer and node 
!        vdens_local     =  ! veg. density / per layer and node 
!        lthick_local    =  ! lay. thickness / per layer and node 
!
        vegdiss(:,IP)   = 0
        vdrgcoeff(:,IP) = 1.
        vdiam(:,IP)     = 0.04
        vdens(:,IP)     = 0.6 
        lthick(:,IP)    = 2.

        CALL INTVEGDISSIP(vegdiss,nvrt,dep(ip),kbar,vdrgcoeff(:,IP),vdiam(:,IP),vdens(:,IP),lthick(:,IP)) 

        SVEGET = BGDISS * VEGDISS 

        DO IS = 1, NUMSIG
          DO ID = 1, NUMDIR
            DSSVEG(IS,ID) = - SVEGET
            SSVEG(IS,ID)  = - SVEGET * ACLOC(IS,ID)
          ENDDO
        ENDDO 
 
      END SUBROUTINE
!****************************************************************
!
!****************************************************************
      subroutine intvegdissip(vegdiss,nlay,depth,kbar,vdrgcoeff,vdiam,vdens,lthick)
        USE DATAPOOL, ONLY : RKIND
        implicit none

        real(rkind),intent(in)  :: depth
        real(rkind),intent(in)  :: kbar
        real(rkind),intent(in)  :: vdrgcoeff(nlay)
        real(rkind),intent(in)  :: vdiam(nlay)
        real(rkind),intent(in)  :: vdens(nlay)
        real(rkind),intent(in)  :: lthick(nlay)
        real(rkind),intent(out) :: vegdiss
        real(rkind)             :: svkh1, svkh2, coeff, kvh

        integer,intent(in)      :: nlay
        integer                 :: i,j

        svkh1 = 0.d0
        svkh2 = 0.d0
        kvh   = 0.d0
        do i = 1, nlay
          sumlay  = sumlay + lthick(i)  
          if (vdiam(i) .gt. ZERO) then
            kvh     = kvh + kbar * lthick(i)
            svkh1   = svkh2 
            svkh2   = svkh2 + sinh(kvh)
            coeff   = (svkh2**3-svkh1**3)+3*(svkh2-svkh1)
            vegdiss = vegdiss + coeff*vdiam(i)*vdens(i)*lthick(i)
          endif
        enddo

      end subroutine
!****************************************************************
!
!****************************************************************
