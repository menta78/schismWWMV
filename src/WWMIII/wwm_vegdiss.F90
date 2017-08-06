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

        INTEGER                 :: IS, ID

        REAL(rkind)             :: BGDISS, KDBAR, SVEGET, VEGDISS
        REAL(rkind)             :: VCD, VDM, VDENS, VLTH(NVRT)

        IF (ETOT .LT. THR) RETURN

        KDBAR  = KBAR * DEPTH
!
        BGDISS = SQRT(TWO/PI)*G9**2*(KBAR/SBAR)**3*SQRT(ETOT)/(THREE*KBAR*(COSH(KDBAR))**3)*VDENS
!
!       VCD =  ! drag coefficient / per layer and node 
!       VDM =  ! diam of veg. / per layer and node 
!       VDENS =  ! veg. density / per layer and node 
!       VLTH =  ! lay. thickness / per layer and node 
!
        CALL INTVEGDISSIP(vegdiss,nvrt,depth,kbar,vcd,vdm,vdens,vlth) 

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
      SUBROUTINE INTVEGDISSIP(vegdiss,nlay,depth,kbar,vdrgcoeff,vdiam,vdens,lthick)
        USE DATAPOOL, ONLY : RKIND, ZERO
        implicit none

        real(rkind),intent(in)  :: depth
        real(rkind),intent(in)  :: kbar
        real(rkind),intent(in)  :: vdrgcoeff(nlay)
        real(rkind),intent(in)  :: vdiam(nlay)
        real(rkind),intent(in)  :: vdens(nlay)
        real(rkind),intent(in)  :: lthick(nlay)
        real(rkind),intent(out) :: vegdiss
        real(rkind)             :: svkh1, svkh2, coeff, kvh, sumlay 

        integer,intent(in)      :: nlay
        integer                 :: i,j

        svkh1 = ZERO
        svkh2 = ZERO
        kvh   = ZERO
        sumlay = ZERO
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
