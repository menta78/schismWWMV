      SUBROUTINE bedchange_bedload(ised,it,moitn,mxitn,rtol,qsan,    &
      &                            bc_sed,lbc_sed,hbed,hbed_ised)    
!--------------------------------------------------------------------!
! This computes changes in bed characteristics and thickness induce  !
! by bedload transport                                               !
!                                                                    !
! Author: Knut Kraemer                                               !
! Date: 27/11/2012                                                   !
!                                                                    !
! History: 2012/12 - F.Ganthy : form homogenisation of sediments     !
!          routines                                                  !
!                                                                    !
!--------------------------------------------------------------------!

      USE sed_mod
      USE elfe_glbl, ONLY : rkind,nm,nea,npa,nne,idry,idry_e,xctr,   &
                            yctr,np,ine,errmsg,xnd,ynd
      USE elfe_msgp
      
      IMPLICIT NONE

!- Local variables --------------------------------------------------!

      LOGICAL,INTENT(IN) :: lbc_sed(npa)  !b.c. flag for erosion eq.,
                                          ! used in JCG
      INTEGER, INTENT(IN) :: ised,it ! sediment class and time step
      INTEGER, INTENT(IN) :: moitn,mxitn ! JCG solver iterations
      REAL(rkind),INTENT(IN) :: rtol     ! Relative tolerance

      INTEGER     :: i,j,nm1,nm2,nm3,ks,k
      REAL(rkind) :: yp,xp,flux,cff,cff1
      !b.c. for erosion eq., used in JCG
      REAL(rkind),DIMENSION(npa) :: bed_poro,qsan,hbed,hbed_ised,    &
                                    bc_sed  
      
!- Start Statement --------------------------------------------------!
      
      WRITE(16,*)'SED: Entering bedchange_bedload'
      
!---------------------------------------------------------------------
! -  qsan is overall sand flux at each node
! qsaxy is the sand flux at the element center integrated in time. Now 
! compute the line integral of qsaxy*normal along the control volume. 
! Add for each node in qsan. The unit normal is directed outward.
!---------------------------------------------------------------------
      qsan=0.0d0

      DO i = 1,nea
        IF(idry_e(i)==1) CYCLE
            
!---------------------------------------------------------------------
! -  Apply morphology factor to bedload transport
!---------------------------------------------------------------------
        FX_r(i) = FX_r(i)*morph_fac(ised)*bed_frac(1,i,ised)
        FY_r(i) = FY_r(i)*morph_fac(ised)*bed_frac(1,i,ised)

        nm1 = nm(i,1)
        nm2 = nm(i,2)
        nm3 = nm(i,3)

!---------------------------------------------------------------------
! - Integrate bedload flux from element center to edge centers
! balance flux at neighbouring edges to obtain flux for nodes
! flux [m3], qsan [m3]
!---------------------------------------------------------------------

        xp   = (xnd(nm2)+xnd(nm3))/2.d0
        yp   = (ynd(nm2)+ynd(nm3))/2.d0
        flux = FX_r(i)*(yp-yctr(i))-FY_r(i)*(xp-xctr(i))
        qsan(nm2) = qsan(nm2)-flux
        qsan(nm3) = qsan(nm3)+flux

        xp   = (xnd(nm3)+xnd(nm1))/2.d0
        yp   = (ynd(nm3)+ynd(nm1))/2.d0
        flux = FX_r(i)*(yp-yctr(i))-FY_r(i)*(xp-xctr(i))
        qsan(nm3) = qsan(nm3)-flux
        qsan(nm1) = qsan(nm1)+flux

        xp   = (xnd(nm1)+xnd(nm2))/2.d0
        yp   = (ynd(nm1)+ynd(nm2))/2.d0
        flux = FX_r(i)*(yp-yctr(i))-FY_r(i)*(xp-xctr(i))
        qsan(nm1) = qsan(nm1)-flux
        qsan(nm2) = qsan(nm2)+flux

      ENDDO !End loop nea

!---------------------------------------------------------------------
! - Compute erosion rates
! change bed(1,mne,iporo) from elements to nodes 
! initalize before adding
!---------------------------------------------------------------------
      bed_poro = 0.0d0

      DO i=1,np
        IF(idry(i)==1) CYCLE

        ks=0 !Number of wet neighbor elements
        DO j=1,nne(i)
          k = ine(i,j)
          IF(idry_e(k)==1) CYCLE
          ks = ks+1
          bed_poro(i) = bed_poro(i)+bed(1,k,iporo)
        ENDDO ! End loop nne
        
        IF(ks==0) CALL parallel_abort('SEDIMENT: (2)')
        bed_poro(i) = bed_poro(i)/ks
      ENDDO ! End loop np

      ! Exchange ghosts
      call exchange_p2d(bed_poro(:))

!...RHS
!      call exchange_p2d(qsan)

!---------------------------------------------------------------------
! - Take porosity into account and adjust qsan accordingly
!jl. mcoefd      [m2] ----  aux1/aux2 are related to Exner equation 
!    hdbed_ised  [m]
!    qsan        [m3]
!    A      x          = B
!    mcoefd hdbed_ised = qsan
!    [m2]   [m]        = [m3]
!---------------------------------------------------------------------
!jl. Why is qsan divided by (1-porosity)?  Doesn't this magically
!   dd volume? I think the idea is to scale the value by (1-porosity).

      qsan(1:np) = qsan(1:np)/(1-bed_poro(1:np))

      ! Use JCG solver to calc hbed_ised
      hbed_ised=0.0d0 !initial guess
      CALL solve_jcg(it,moitn,mxitn,rtol,mcoefd,hbed_ised,qsan,      &
      &              bc_sed,lbc_sed)

!---------------------------------------------------------------------
! - Bed/bottom level change due to bedload transport in [m]
! ????? isn't there a dimension mismatch ?????
!---------------------------------------------------------------------

      hbed(:) = hbed(:)+hbed_ised(:)


      ! Consistency check
      DO i=1,np
        IF (isnan(hbed(i)).EQV..TRUE.) THEN
          WRITE(errmsg,*)'hbed(i) is NaN',myrank,i,hbed(i),qsan(i),  &
          &              bed_poro(i)
          CALL parallel_abort(errmsg)
        ENDIF
      ENDDO

!---------------------------------------------------------------------
! - Evaluate depth at the elements and changes in bed properties
!---------------------------------------------------------------------
      DO i=1, nea
        IF (idry_e(i)==1) CYCLE

        ! Average bed change in element [m]
        cff = (hbed_ised(nm(i,1))+                                   &
        &      hbed_ised(nm(i,2))+                                   &
        &      hbed_ised(nm(i,3)))/3.d0

        ! Average bed change in element [kg/m2]
        cff1=cff*Srho(ised)

!---------------------------------------------------------------------
! - Update bed mass [kg/m2] according to bed change
!jl. You can lose a maximum of bed_mass per sediment class per time 
!    step based on changes to hbed_ised.  Is there a lower limit on 
!    hbed_ised?...
!    No lower limit on hbed_ised, but bed_mass and bed(1, nea, ithick)     
!    have a lower bound of 0.
!---------------------------------------------------------------------

        bed_mass(1,i,nnew,ised) = MAX(bed_mass(1,i,nstp,ised) +      &
        &                             cff1,0.0d0)

       !Jan what is this for?
        IF (suspended_load == 1) THEN
          DO k=2,Nbed
            bed_mass(k,i,nnew,ised) = bed_mass(k,i,nstp,ised)
          ENDDO
        ENDIF

!---------------------------------------------------------------------
! - Update layer thickness according to bed change in [m]
!---------------------------------------------------------------------

        bed(1,i,ithck) = MAX((bed(1,i,ithck)+cff),0.0d0)
      ENDDO !End loop nea


!--------------------------------------------------------------------!
      END SUBROUTINE bedchange_bedload
