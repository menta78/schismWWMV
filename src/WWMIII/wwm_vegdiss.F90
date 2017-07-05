!****************************************************************
!
!****************************************************************
#ifdef SCHISM
      SUBROUTINE VEGDISSIP (SDVEG,DDVEG,ACLOC,DEPTH,ETOT,SBAR,KBAR) 
        USE DATAPOOL
        IMPLICIT NONE
        REAL(rkind),INTENT(IN)  :: KBAR,SBAR,ETOT,DEPTH
        REAL(rkind),INTENT(IN)  :: ACLOC(NUMSIG,NUMDIR)
        REAL(rkind),INTENT(OUT) :: SDVEG(NUMSIG,NUMDIR), DDVEG(NUMSIG,NUMDIR)

        REAL(rkind) :: BGDISS, KDBAR, SVEGET, VEGDISS
        REAL(rkind) :: vdrgcoeff_local(nvrt), vdiam_local(nvrt), vdens_local(nvrt), lthick_local(nvrt)
        INTEGER :: IS, ID

        KDBAR  = KBAR * DEPTH
!
        BGDISS = SQRT(TWO/PI)*G9**2*(KBAR/SBAR)**3*SQRT(ETOT)/(THREE*KBAR*(COSH(KDBAR))**3)!*NPLANTSPSQM(IP)
!
!2do ... extract from schism 

!        vdrgcoeff_local =  ! drag coefficient / per layer and node 
!        vdiam_local     =  ! diam of veg. / per layer and node 
!        vdens_local     =  ! veg. density / per layer and node 
!        lthick_local    =  ! lay. thickness / per layer and node 
!
        IF (LVEGCONST) THEN
#ifdef SCHISM
          CALL intvegdissip_const_schism(vegdiss,vdrgcoeff_local,vdiam_local,vdens_local,kbar,lthick_local)
#else
          CALL intvegdissip_const(vegdiss,vdrgcoeff_local,vdiam_local,vdens_local,kbar,lthick_local) 
#endif
        ELSE
          CALL intvegdissip_var(vegdiss,nvrt,vdrgcoeff_local,vdiam_local,vdens_local,kbar,lthick_local) 
        ENDIF
        
        SVEGET = BGDISS + VEGDISS 
    
        DO IS = 1, NUMSIG
          DO ID = 1, NUMDIR
            DDVEG(IS,ID) = - SVEGET
            SDVEG(IS,ID) = - SVEGET * ACLOC(IS,ID)
          ENDDO
        ENDDO 

      END SUBROUTINE
!****************************************************************
!
!****************************************************************
      subroutine intvegdissip(vegdiss,nlay,vdrgcoeff,vdiam,vdens,kbar,lthick)
        USE DATAPOOL, ONLY : RKIND
        implicit none
        real(rkind),intent(in)  :: kbar
        real(rkind),intent(in)  :: vdrgcoeff(nlay)
        real(rkind),intent(in)  :: vdiam(nlay)
        real(rkind),intent(in)  :: vdens(nlay)
        real(rkind),intent(in)  :: lthick(nlay)
        real(rkind),intent(out) :: vegdiss
        integer,intent(in)      :: nlay
        integer                 :: i,j
        real(rkind)             :: svkh1, svkh2, sum_thick, coeff
        svkh1 = 0.d0
        svkh2 = 0.d0
        sum_thick = 0.d0
        do i = 1, nlay
          svkh1            = svkh1 + sinh(kbar * lthick(i-1))
          if(i.gt.1) svkh2 = svkh2 + sinh(kbar * lthick(i))
          coeff            = (svkh2**3-svkh1**3)+3*(svkh2-svkh1)
          vegdiss          = vegdiss + coeff*vdiam(i)*vdens(i)*lthick(i)
        enddo
      end subroutine
#endif
!****************************************************************
!
!****************************************************************
