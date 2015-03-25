      SUBROUTINE sed_roughness
!--------------------------------------------------------------------!
! This routine computes ripples and sand-waves parameters, and then  !
! update the roughness for hydrodynamics (rough_p) length and for    !
! sediment model (Zob) accordingly.                                  !
! By setting bedforms_rough to 1, 2 or 3 within sediment.in, you can !
! choose to update both hydro and sediment, only hydro or only       !
! sediment related roughness length                                  !
!                                                                    !
! Author: florian ganthy (fganthy@lnec.pt ; florian.ganthy@gmail.com)!
! Date:   2013/01/02                                                 !
!                                                                    !
! History:                                                           !
!                                                                    !
!--------------------------------------------------------------------!
      USE elfe_glbl, ONLY: npa,rho0,pi,grav,kbp,znl,uu2,vv2,idry,h0, &
                           eta2,dpe,errmsg,nm,nea,rkind,rough_p,dfv
      USE elfe_msgp, ONLY: myrank,parallel_abort
      USE sed_mod,   ONLY: bed_d50n,bed_taucen,bed_ripplen,          &
                           bed_ripphgt,bed_rhosn,bed_z0defn,Zob,     &
                           vonKar,rough_ripple,rough_sedtrans,       &
                           Cdb_min,Cdb_max,bedforms_rough,bed_rough

      IMPLICIT NONE
      SAVE

!- Local variables --------------------------------------------------!

      INTEGER     :: i,n1,n2,n3
      REAL(rkind) :: hwat,z,cc1,cc2,cc3,wrk
      REAL(rkind) :: tau_c,tau_cr,tau_wash,tau_s
      REAL(rkind) :: z0nik,z0bdf,z0st
      REAL(rkind) :: ar,a0,a1,a2

!- Start Statement --------------------------------------------------!

      IF(myrank.EQ.0) WRITE(16,*)'Entering sed_roughness'

      bed_rough = bed_z0defn
      ! coefficient for bedforms roughness (Nielsen, 1992)
      ar    = 0.267d0

      DO i=1,npa

        hwat   = dpe(i)+eta2(i) ! Total water depth (m)

!---------------------------------------------------------------------
! ** Nikuradse roughness (independant from wet-drying)
!---------------------------------------------------------------------
          z0nik = bed_d50n(i)/12.d0 ! in (m)


        IF(idry(i).EQ.1.OR.hwat.LE.h0) THEN
          ! Dry nod
          z0st = 0.d0 !Sediment transport roughness length
  
        ELSE
          ! Wet nod

          ! Total Critical shear stress (depends on sediment
          ! composition)
          tau_cr = bed_taucen(i)  ! tau_cr in m2.s-1

!---------------------------------------------------------------------
! ** Current-induced bed shear stress (skin friction only)
!---------------------------------------------------------------------
          z = znl(kbp(i)+1,i)-znl(kbp(i),i) ! Height above bed (m)
          IF(z.GT.z0nik) THEN
            ! Bed shear stress for logarithmic part of the BBL
            IF((z.LE.0).OR.(z0nik.LE.0)) THEN
              WRITE(errmsg,*)'SED ROUGHNESS: Cd Failed,',myrank,i,z, &
              &              z0nik
              CALL parallel_abort(errmsg)
            ENDIF
            cc1   = 1.d0/DLOG(z/z0nik)
            cc2   = vonKar**2.d0*cc1**2.d0
            wrk   = MIN(Cdb_max,MAX(Cdb_min,cc2))
            tau_c = wrk*uu2(kbp(i)+1,i)**2.d0*vv2(kbp(i)+1,i)**2.d0
          ELSE
            ! Bed shear stress for the viscous sublayer
            z = znl(kbp(i)+2,i)-znl(kbp(i)+1,i)
            IF(z.LE.0) THEN
              WRITE(errmsg,*)'SED ROUGHNESS: div. by 0,',myrank,i,z, &
              &              z0nik
              CALL parallel_abort(errmsg)
            ENDIF
            cc1   = (uu2(kbp(i)+2,i)-uu2(kbp(i)+1,i))
            cc2   = (vv2(kbp(i)+2,i)-vv2(kbp(i)+1,i))
            cc3   = (dfv(i,kbp(i)+1)-dfv(i,kbp(i)+2))/2
            tau_c = SQRT((cc1*cc3/z)**2.d0+(cc2*cc3/z)**2.d0)
          ENDIF

!---------------------------------------------------------------------
! ** Current-induced bedforms
!---------------------------------------------------------------------
          IF(tau_c.GT.tau_cr)THEN

            ! Enough bed shear stress for sediment transport
            IF(rough_ripple.EQ.1) THEN

              ! Ripples (Soulsby, 1999)
              tau_wash = 0.8d0*grav*((bed_rhosn(i)/rho0) -1.d0)*     &
              &          bed_d50n(i)
              IF(tau_c.LT.tau_wash) THEN
                bed_ripplen(i) = 1000.d0*bed_d50n(i)     ! Rip. length
                bed_ripphgt(i) = bed_ripplen(i)/7.d0     ! Rip. height
              ELSE
                bed_ripphgt(i) = 0.d0        ! Rip. height
                bed_ripplen(i) = 1.d0        ! Rip. length (arbitrary)
              ENDIF ! Test tau_c<tau_wash

            ELSEIF(rough_ripple.EQ.2) THEN

              ! Sand waves (van Rijn, 1984):
              tau_wash = tau_cr*26.d0
              IF(tau_c.LT.tau_wash) THEN
                tau_s          = (tau_c-tau_cr)/tau_cr
                bed_ripphgt(i) = 0.11d0*hwat*                        &
                &                ((bed_d50n(i)/hwat)**0.3d0)*        &
                &                (1.d0-EXP(1.d0)**(-0.5d0*tau_s))*   &
                &                (25.d0-tau_s)         ! Rip. height
                bed_ripplen(i) = 7.3d0*hwat            ! Rip. lenght
              ELSE
                bed_ripphgt(i) = 0.d0       ! Rip. height
                bed_ripplen(i) = 1.d0       ! Rip. length (arbitrary)
              ENDIF ! Test tau_c<tau_wash

            ELSEIF(rough_ripple.EQ.3) THEN

              ! Sand waves (Yalin, 1964)
              tau_wash = tau_cr*17.6d0
              IF(tau_c.LT.tau_wash) THEN
                tau_s          = tau_c/tau_cr
                bed_ripphgt(i) = (hwat/6.d0)*(1.d0-tau_s) !Rip. height
                bed_ripplen(i) = 2.d0*pi*hwat             !Rip. length
              ELSE
                bed_ripphgt(i) = 0.d0        ! Rip. height
                bed_ripplen(i) = 1.d0        ! Rip. length (arbitrary)
              ENDIF

            ELSE

              ! No taking account for current ripples
              bed_ripphgt(i) = 0.d0          ! Rip. height
              bed_ripplen(i) = 1.d0          ! Rip. length (arbitrary)
            ENDIF ! END rough_ripple (Soulsby or van Rijn)
          ELSE

            ! No sediment transport, bedforms from previous time step
            bed_ripphgt(i) = bed_ripphgt(i)              ! Rip. Height
            bed_ripplen(i) = bed_ripplen(i)              ! Rip. length
          ENDIF ! END tau_c>tau_cr


          IF(isnan(bed_ripphgt(i)).EQV..TRUE.) THEN
            WRITE(errmsg,*)'SED ROUGHNESS: ripple height NaN',myrank,&
           &                i,bed_ripphgt(i),tau_c,bed_d50n(i),hwat
           CALL parallel_abort(errmsg)
          ENDIF
          IF(bed_ripphgt(i).LT.0.d0) THEN
            WRITE(errmsg,*)'SED ROUGHNESS: ripple height <0',myrank,i, &
            &                bed_ripphgt(i),tau_c,bed_d50n(i),hwat
            CALL parallel_abort(errmsg)
          ENDIF         

!---------------------------------------------------------------------
! ** Sediment transport roughness (Warner et al., 2008):
!---------------------------------------------------------------------
          ! Coefficients (Warner et al., 2008, Wiberg and Rubin, 1989)
          a0   = 0.056d0
          a1   = 0.068d0
          a2   = 0.0204d0*LOG(100.d0*(bed_d50n(i)**2.d0))+             &
          &      0.0709d0*LOG(100.d0*bed_d50n(i))
          tau_s = tau_c/tau_cr

          IF(rough_sedtrans.EQ.1.AND.tau_s.GT.1.)THEN
            z0st = a0*bed_d50n(i)*a1*tau_s/(1.d0+a2*tau_s)
          ELSE
            z0st = 0.d0
          ENDIF


        ENDIF ! End test on wet-dry nod

!---------------------------------------------------------------------
! ** Bedforms roughness:
!---------------------------------------------------------------------

        IF(rough_ripple.EQ.0)THEN
          z0bdf = 0.d0
        ELSE
          z0bdf = ar*(bed_ripphgt(i)**2.d0)/bed_ripplen(i)
        ENDIF

!---------------------------------------------------------------------
! ** Total roughness length
!---------------------------------------------------------------------

        bed_rough(i) = MAX(z0nik+z0bdf+z0st,bed_z0defn(i))

        IF(bed_rough(i).LE.0.d0) THEN
         WRITE(errmsg,*)'SED ROUGHNESS: bed_rough(i) <= 0',myrank,i, &
         &              bed_rough(i)
         CALL parallel_abort(errmsg)
        ENDIF
        IF(bed_rough(i).GE.0.5d0) THEN
          WRITE(errmsg,*)'SED ROUGHNESS: To high roughness! >0.5m ',  &
          myrank,i,bed_rough(i),z0nik,z0bdf,z0st,bed_ripphgt(i),      &
          bed_ripplen(i),bed_d50n(i)
          CALL parallel_abort(errmsg)
        ENDIF

!---------------------------------------------------------------------
! ** Applying total roughness length to hydrodynamics roughness
!---------------------------------------------------------------------

        IF(bedforms_rough.EQ.1.OR.bedforms_rough.EQ.2.OR.bedforms_rough.EQ.4)THEN
          rough_p(i) = bed_rough(i)
        ENDIF

      ENDDO !END loop npa

!---------------------------------------------------------------------
! ** Applying total roughness length to element for sediment transport
!---------------------------------------------------------------------

      IF(bedforms_rough.EQ.1.OR.bedforms_rough.EQ.3)THEN
        ! Application of total bedforms roughness length
        DO i=1,nea
          n1     = nm(i,1)
          n2     = nm(i,2)
          n3     = nm(i,3)
          Zob(i) = MAX((bed_rough(n1)+                               &
          &             bed_rough(n2)+                               &
          &             bed_rough(n3))/3.d0,1.d-5)
          IF(Zob(i).LE.0.d0) THEN
             WRITE(errmsg,*)'SED ROUGHNESS: Zob(i) <= 0',myrank,i,   &
             &              Zob(i)
             CALL parallel_abort(errmsg)
          ENDIF
        ENDDO
      ELSEIF(bedforms_rough.EQ.4)THEN
        ! Application of nikuradse roughness length
        DO i=1,nea
          n1     = nm(i,1)
          n2     = nm(i,2)
          n3     = nm(i,3)
          Zob(i) = MAX((bed_d50n(n1)+                                &
          &             bed_d50n(n2)+                                &
          &             bed_d50n(n3))/36.d0,1.d-5)
          IF(Zob(i).LE.0.d0) THEN
             WRITE(errmsg,*)'SED ROUGHNESS: Zob(i) <= 0',myrank,i,   &
             &              Zob(i)
             CALL parallel_abort(errmsg)
          ENDIF
        ENDDO
      ENDIF

!---------------------------------------------------------------------
      IF(myrank.EQ.0) WRITE(16,*)'Leaving sed_roughness'
!--------------------------------------------------------------------!
      END SUBROUTINE sed_roughness   
