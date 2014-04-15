!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine init_flmud

        use elfe_glbl, only : npa, nvrt, kbp, rkind, irouse_test,idry
        use flmud_pool
        implicit none

        integer :: il,ip
!        integer :: testnode = 1000

        real(rkind) :: phi(nvrt), phip(nvrt), phigel(nvrt), cgel(nvrt) !, nf(nvrt)
!        logical       :: ldebug = .true.

!       Read in parameters
        open(31,file='timor.in',status='old')
        read(31,*)
        read(31,*)ithick
        read(31,*)dm0
        read(31,*)dm1
        read(31,*)mudref
        read(31,*)laddmud_d
        read(31,*)laddmud_v
        read(31,*)ldumpmud
        read(31,*)lhindws
        read(31,*)lfloc
        read(31,*)ldebug
        read(31,*)testnode
!        read(31,*)icycle
        close(31)

!       --------------------------------------------------
!       allocate mud stuff ... 
!       --------------------------------------------------
        call init_flmud_arrays
!        write(12,*)'done init mud arrays'
!       --------------------------------------------------
!       init mud stuff ... 
!       --------------------------------------------------
        dmf_mud = dm1 !100./1.e6 !10**6.
        nfg=2 !init. for dry nodes as well
        do ip = 1, npa

          if(idry(ip)==1) cycle
!           --------------------------------------------------
!           add mud at first time step 
!           --------------------------------------------------
!          nfg(:,ip) = 2.
          do il = kbp(ip),nvrt
            if (il-kbp(ip)+1<=ithick) then
              mudc(1,il,ip) = mudref
            elseif(il-kbp(ip)+1>ithick) then
              mudc(1,il,ip) = 0
            end if
          end do
!          write(12,*)'done init mudc ',ip

!         --------------------------------------------------
!         Set mud density 
!         --------------------------------------------------
          call get_mudrho(ip,nfg(:,ip))
          !write(12,*)'done init get_mudrho',ip
!         --------------------------------------------------
!         Set gelling concentration
!         --------------------------------------------------
          call get_cgel(ip,nfg(:,ip),cgel,phigel)
!          write(12,*)'done init get_cgel',ip
!         --------------------------------------------------
!         Set mud vol. concentration
!         --------------------------------------------------
          call get_mudc(ip,nfg(:,ip),phi,phip)
!          write(12,*)'done init get_mudc',ip
!         --------------------------------------------------
!         Set sink vel. ... this is overwritten ... 
!         --------------------------------------------------
          call set_hind_wsink(ip,phi,phip,phigel)
!          write(12,*)'done init set_hind_wsink',ip
          !wsink(:,k) = 0.
       
!         Rouse test
          if(irouse_test==1) then
            wsink(:,:,ip)=5.e-2
            !The following actually not necessary as trel will be overwritten
            mudc(:,:,ip)=0
            mudc(:,kbp(ip)+1,ip)=1
            mudc(:,kbp(ip),ip)=1
          endif !irouse_test
!          write(12,*)'done node:',ip
        end do !ip
!        write(12,*)'done init mud...'

      end subroutine init_flmud

!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine flmud(ip,dt,rough_p,SS1d,NN1d,tke1d,eps1d,L1d,num1d,nuh1d)

        use elfe_glbl, only : npa, nvrt, kbp, rkind, irouse_test
        use flmud_pool 
        implicit none

        integer :: il

        integer, intent(in) :: ip
        real(rkind), intent(in) :: dt
        real(rkind), intent(in) :: rough_p !(npa)
        real(rkind), intent(in), dimension(0:nvrt) :: SS1d,NN1d,tke1d,eps1d,L1d,num1d,nuh1d


        real(rkind) :: phi(nvrt), phip(nvrt), phigel(nvrt), cgel(nvrt) !, nf(nvrt)
        real(rkind) :: fak, smooth

        logical, save :: linit
        real(rkind), save :: t_start, t_now

        data t_now /1./
        data linit /.true./

        t_now = t_now + 1.
        t_start = 1000.
        smooth  = 100.

!       Ramp up mud conc.
        if (ldumpmud) then
!          do ip = 1,npa
!         --------------------------------------------------
!         add mud linear with time at a certain t_start time
!         --------------------------------------------------
          do il = kbp(ip), nvrt
            if (il-kbp(ip)+1>ithick) then
              mudc(1,il,ip) = 0.
            elseif(il-kbp(ip)+1<=ithick) then
              fak = 0.5d0+0.5d0*tanh((t_now-t_start)/smooth)
              if (fak .gt. 0.9999999) cycle
              mudc(1,il,ip) = fak * mudref + (1.-fak)*0.
              !if(k==103)write(*,*) fak, mudc(l,k)
            end if
          end do
!          end do !ip
        endif !ldumpmud

!       Time loop
!        do ip = 1, npa
!       ------------------------------------------------------
!       compute new floc diameter
!       ------------------------------------------------------
        if (lfloc) then
          call get_floc_diam(ip,dt,nfg(:,ip))
        else
          nfg(:,ip) = 2.
        endif
!       ------------------------------------------------------
!       Get mud density  
!       -----------------------------------------------------
        call get_mudrho(ip,nfg(:,ip))
!       ------------------------------------------------------
!       Get gelling concentration  
!       -----------------------------------------------------
        call get_cgel(ip,nfg(:,ip),cgel,phigel)
!       ------------------------------------------------------
!       Get mud vol. concentration  
!       ------------------------------------------------------
        call get_mudc(ip,nfg(:,ip),phi,phip)
!       -------------------------------------------------------
!       Update settling velocity
!       -------------------------------------------------------
        call set_hind_wsink(ip,phi,phip,phigel)
!       -------------------------------------------------------
!       Update rheological stress and viscosity
!       -------------------------------------------------------
        call stress_mud(ip,SS1d)
!       -------------------------------------------------------
!       Write test output if asked so ... 
!       -------------------------------------------------------
        call write_point_output(ip,SS1d,NN1d,tke1d,eps1d,L1d,num1d,nuh1d)
 
        end subroutine flmud
!**********************************************************************
!*                                                                    *
!**********************************************************************
        subroutine init_flmud_arrays
          
          use elfe_glbl, only: nvrt,npa,ntracers,rhosed
          use flmud_pool
           
          implicit none

          allocate( rhomud(ntracers,nvrt,npa) ); rhomud = rhosed !for init. call of eqstate and dry nodes
          allocate( dmf_mud(ntracers,nvrt,npa) ); dmf_mud = dm1
          allocate( nfg(nvrt,npa) )
          allocate( wsink(ntracers,nvrt,npa) ); wsink = 0.
          allocate( vts(nvrt,npa) ); vts = 0.
          allocate( dfv_yield(nvrt,npa) ); dfv_yield = 0.
          allocate( dfh_yield(nvrt,npa) ); dfh_yield = 0.
          allocate( rough_mud(npa) ); rough_mud = 0.

       end subroutine
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine get_mudc(ip,nf,phi,phip)
        use elfe_glbl, only: kbp, nvrt, rkind,iplg
        use flmud_pool 

        implicit none

        integer, intent(in)            :: ip
        real(rkind), intent(in)        :: nf(nvrt)      ! Fractal dimension
!        integer, intent(in)            :: testnode
        real(rkind), intent(out)       :: phip(nvrt)
        real(rkind), intent(out)       :: phi(nvrt)
!        logical, intent(in)            :: ldebug
        integer                        :: il   

        do il = kbp(ip), nvrt
          phip(il) = mudc(1,il,ip) /rhomud(1,il,ip) !rhosed
          !phi(il)  = mudc(1,il,ip) / rhosed * (dmf_mud(1,il,ip)/dm0)**(3-nf(il)) ! nvrt npa
          phi(il)  = mudc(1,il,ip)/rhomud(1,il,ip)* (dmf_mud(1,il,ip)/dm0)**(3-nf(il)) ! nvrt npa
          !phi(il)=((rhosed-rhow)/(rhomud(1,il,ip)-rhow))*(mudc(1,ip,il)/rhosed)
          if (ldebug .and. iplg(ip) == testnode) then
            write(1111,'(A10,I10,6F16.8)') 'get_mudc', il,              &
     &                                phip(il), phi(il),                &
     &                                rhosed, rhow,                     &
     &                                rhomud(:,il,ip), mudc(1,il,ip)
          endif
        end do

       end subroutine get_mudc
!**********************************************************************
!*                                                                    *
!**********************************************************************
       subroutine get_mudrho(ip,nf)

        use elfe_glbl, only: idry, kbp, nvrt, rkind,iplg
        use flmud_pool 
        implicit none


        integer, intent(in)                 :: ip
!        integer, intent(in)                 :: testnode
        integer                             :: il
        real(rkind), intent(in)        :: nf(nvrt)         ! Fractal dimension
!        logical, intent(in)                 :: ldebug

        do il = KBP(ip), NVRT
          rhomud(1,il,ip) = rhow + (rhosed - rhow) * (dm0/dmf_mud(1,il,ip))**(3.-nf(il))
           if (ldebug .and. iplg(ip) == testnode) then
             write(1112,'(A10,I10,7F15.8)') 'get_mudrho', il,              & 
     &                  rhow, rhosed, rhomud(1,il,ip),                       &
     &                  dm0, dmf_mud(1,il,ip),nf(il)
           endif
        enddo

       end subroutine get_mudrho
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine set_hind_wsink(ip,phi,phip,phigel)
!
        use elfe_glbl, only: nvrt, kbp, grav, rkind,iplg
        use flmud_pool
!
        implicit none

        real(rkind)                    :: d              !grain size
        real(rkind)                    :: visk           !molecular viskosity 
        integer                        :: il
        real(rkind)                    :: rhost
        real(rkind)                    :: nu
        real(rkind)                    :: dst
        real(rkind)                    :: ws0
        real(rkind)                    :: ws
!        real(rkind)                    :: rhow
        real(rkind)                    :: fak
        real(rkind)                    :: dmi
        real(rkind)                    :: tf

        integer, intent(in)            :: ip
        real(rkind), intent(in)        :: phi(nvrt)
        real(rkind), intent(in)        :: phip(nvrt)
        real(rkind), intent(in)        :: phigel(nvrt)
!        integer, intent(in)            :: testnode
!        logical, intent(in)            :: ldebug

        do il = kbp(ip), nvrt
           !rhost = (rhosed-rhow)/rhow ! 1.65
           rhost = (rhomud(1,il,ip)-rhow)/rhow ! 1.65
           dmi   = dmf_mud(1,il,ip)
           nu    = 1.73E-6 !vismol!max(vismol,dfv(ip,il)) ! 1.E-6
           dst   = (rhost*grav/nu**2)**(1./3.)*dmi ! approx. 3.5 for particles of 200micros
           ws0   = 11*nu/dmi*(SQRT(1+.01*dst**3)-1)
           fak   = (1.-min(1.d0,phi(il))*(1.-phip(il)))/(1+2.5*phi(il))
           if(lhindws) then
             ws=ws0*fak
           else
             ws=ws0 !Stoke's law only
           endif
!           wsink(1,il,ip) = ws 
!           if(phi(il).gt.phigel(il)) then
!             wsink(1,min(nvrt,il+1),ip) = 0.
!           tf = 0.5d0+0.5d0*dtanh((phi(il)-phigel(il))/1.e-12)
!           wsink(1,min(nvrt,il+1),ip)= (1.-tf)*ws
!           endif
           if(il==kbp(ip)) then 
             wsink(1,il,ip) = 0.
             ws             = 0.
!nvrt-1,nvrt may have jump
           else if(il==nvrt) then
             wsink(1,il,ip) = ws
           else
             tf = 0.5d0+0.5d0*dtanh((phi(il)-phigel(il))/1.e-18)
             wsink(1,il,ip)= (1.-tf)*ws
           endif

           if (ldebug .and. iplg(ip)==testnode)then
             write(1114,'(A10,I10,50F16.8)') 'set_wsink',     &
     &         il,                                            &
     &         dmi,                                           &
     &         dst,                                           &
     &         ws0,                                           &
     &         ws,                                            &
     &         fak,                                           &
     &         phi(il),                                       &
     &         phip(il),                                      &
     &         phigel(il), &
     &         wsink(1,il,ip)
             endif
         end do !il

        end subroutine set_hind_wsink
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine stress_mud(ip,SS1d)
        use elfe_glbl, only : grav, kbp, nvrt, rkind, prho,iplg
        use flmud_pool 
        implicit none

        integer, intent(in)                 :: ip
        real(rkind), intent(in)             :: SS1d(0:nvrt)

        integer     :: il

        real(rkind) :: g_dot,g_dot_thr
        real(rkind) :: rhop
        real(rkind) :: visk,tau,viskmax
        real(rkind) :: lambda_e
        real(rkind) :: tf, smooth
        real(rkind) :: mu8,mu0,beta,c,tau0

        g_dot_thr = 0.00001
        smooth    = 1.e-2
        viskmax   = 200.

        do il = kbp(ip), nvrt

           rhop   = max(0.d0,prho(il,ip)/rhow-1.)
           g_dot  = sqrt(SS1d(il-kbp(ip))) !SS1d index translated due to GOTM

           call set_toorman_constants(rhop,mu8,mu0,beta,c)
           call set_yieldstress(rhop,tau0)

           if (g_dot>g_dot_thr.and.rhop>0) then 
             lambda_e = 1.d0/(1.d0+beta*g_dot)
             tau  = tau0 + (mu8 + c + beta*tau0*lambda_e) * g_dot
             visk=-beta**2*tau0*g_dot/(1.d0+beta*g_dot)**2+mu8+c+beta*tau0/(1.+beta*g_dot)
             visk = visk/1000.
           else
             lambda_e  = 0.
             tau       = 0.
             visk      = 0. 
           endif

!           tf = 0.5d0+0.5d0*dtanh((g_dot-g_dot_thr)/smooth)
!           visk = tf * visk + (1.-tf)*viskmax
           vts(il,ip)=visk
 
           if(ldebug .and. iplg(ip) == testnode) then
             write(1115,'(A10,I5,14F15.8)') 'mudvisc',il,rhop,              &
     &                prho(il,ip),g_dot, & !,tf,                                   &
     &                mu0,mu8,beta,c,lambda_e,                              &
     &                tau0,tau,vts(il,ip)
           endif
        end do !il

      end subroutine stress_mud
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine set_yieldstress(rhop,tau_y)
        use elfe_glbl, only: rkind
        implicit none
!
        real(rkind), INTENT(IN)  :: rhop
        real(rkind), INTENT(OUT) :: tau_y

        tau_y = 5000.d0*rhop**3-340.d0*rhop**2+10.d0*rhop

      end subroutine set_yieldstress
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine set_toorman_constants(rhop,mu8,mu0,beta,c)
        use elfe_glbl, only: rkind 
        implicit none
!
        real(rkind), INTENT(IN)  :: rhop
        real(rkind), INTENT(OUT) :: mu8,mu0,beta,c 

        real(rkind)  :: rhop2

        rhop2= rhop*rhop

        mu8  = 1.16*rhop2+0.1*rhop
        mu0  = 0.012*exp(63*rhop)-0.012
        beta = 72*rhop2+6*rhop
        c    = mu0-mu8
        !write(*,'(5F20.10)') rhop, mu8, 1.16*rhop**2.+0.1*rhop

      end subroutine set_toorman_constants
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine set_mud_roughness(ip,alpha,SS1d)
        use elfe_glbl, only: uu2, vv2, kbp, znl, dfv, rkind, prho,nvrt
        use flmud_pool
        implicit none

        integer, intent(in)                 :: ip
        real(rkind), intent(out)             :: alpha
        real(rkind), intent(in)             :: SS1d(0:nvrt)

        real(rkind)                         :: aux,dh,du,dv,m2,dbuoy
        real(rkind)                         :: rho1, rho2, visk_bar,tstress,ufric
        real(rkind)                         :: cnpar,stress_x,stress_y, rhobar
        real(rkind)                         :: g_dot, rhop, tau,visk, drho, ri

        real(rkind), parameter                     :: bpar = 1.
        real(rkind), parameter                     :: beta = 0.7
        real(rkind), parameter                     :: mpar = 1.

        rho1 = prho(kbp(ip),ip) ! is this at the bottom ... selfe starts to count from the bottom?
        rho2 = prho(kbp(ip)+1,ip)
        rhobar = 0.5 * (rho1 + rho2)
        dh = znl(kbp(ip)+1,ip)-znl(kbp(ip),ip)
        drho = abs((rho2-rho1))/dh
!        du = abs((uu2(kbp(ip)+1,ip)-uu2(kbp(ip),ip)))/dh
!        dv = abs((vv2(kbp(ip)+1,ip)-vv2(kbp(ip),ip)))/dh
        g_dot =sqrt(SS1d(0)) !sqrt(du**2+dv**2)
        if (g_dot .gt. 0.0000001) then
          visk_bar = dfv(kbp(ip),ip) !+ vts(kbp(ip),ip)
          stress_x = visk_bar * du
          stress_y = visk_bar * dv
          tstress  = (stress_x**2.+stress_y**2.)**0.5
          ufric = sqrt(tstress)
          ri = 9.81/rhobar*drho/g_dot**2
          alpha = exp(-(1+beta*wsink(1,kbp(ip),ip)/ufric)*(1-exp(-bpar*ri**mpar)))
        else 
          alpha = 1 
        endif
!DS: was ist hier los? alpha überschrieben
          alpha = 1 !

        !write(*,*) ip, alpha, wsink(1,kbp(ip),ip), ri, g_dot**2, ufric

      end subroutine set_mud_roughness
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine set_bottom_stress(ip,bstress,SS1d)
        use elfe_glbl, only: nvrt, kbp, uu2, vv2, dfv, rkind, znl,nvrt
        use flmud_pool
        implicit none

        integer, intent(in)                 :: ip
        real(rkind), intent(out)            :: bstress
        real(rkind), intent(in)            :: SS1d(0:nvrt)

        real(rkind) g_dot, rhop, visk, du, dv, dh, visk_bar
        real(rkind) stress_x, stress_y

!        dh = znl(kbp(ip)+1,ip)-znl(kbp(ip),ip)
!        du = abs((uu2(kbp(ip)+1,ip)-uu2(kbp(ip),ip)))/dh
!        dv = abs((vv2(kbp(ip)+1,ip)-vv2(kbp(ip),ip)))/dh
        g_dot =sqrt(SS1d(0)) !sqrt(du**2+dv**2)
        visk_bar = dfv(kbp(ip),ip) !+ vts(kbp(ip),ip)  ! viscosity
        stress_x = visk_bar * du
        stress_y = visk_bar * dv
        bstress  = (stress_x**2.+stress_y**2.)**0.5
        !write(*,*) ip, visk_bar, tstress, du, dv, dh,kbp(ip)

      end subroutine set_bottom_stress
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine set_bottom_visk(ip,ufric,visk1,visk2)

        use elfe_glbl, only: nvrt, kbp, rkind, znl
        implicit none

        integer                             :: nlev
        integer, intent(in)                 :: ip
        real (rkind), intent(out)           :: visk1
        real (rkind), intent(out)           :: visk2

        real (rkind)                        :: ufric
        real (rkind)                        :: depth
        real (rkind)                        :: dh,h1,h2

        h1 = znl(kbp(ip)+1,ip)-znl(kbp(ip),ip)
        h2 = znl(kbp(ip)+2,ip)-znl(kbp(ip),ip)

        depth  = znl(nvrt,ip)-znl(kbp(ip),ip)
!mixing length ...
        visk1 = 0.41 * ufric * h1*(1 - h2/depth)
        visk2 = 0.41 * ufric * h2*(1 - h2/depth)
!        write(*,'(I10,5F15.10)')ip,visk1,visk2,h(nlev-1),depth,ufric

      end subroutine set_bottom_visk
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine get_floc_diam(ip,dt,nf)
        use elfe_glbl, only: dfv, kbp, nvrt, znl, rkind, errmsg,iplg
        use elfe_msgp, only: parallel_abort
        use flmud_pool
        implicit none

        integer, intent(in)        :: ip
!        integer, intent(in)        :: testnode
!        logical, intent(in)        :: ldebug
        real(rkind), intent(in)    :: dt ! time step
        real(rkind), intent(inout) :: nf(nvrt) ! 3D fractal dimension of floc

        real(rkind)                :: dddt
        real(rkind)                :: dddt1
        real(rkind)                :: dddt2
        real(rkind)                :: dddt3
        real(rkind)                :: beta
        real(rkind)                :: dv
        real(rkind)                :: du
        real(rkind)                :: dh
        real(rkind)                :: cnpar
        real(rkind)                :: g_dot_thr
        real(rkind)                :: g_dot
        real(rkind)                :: g2_dot
        real(rkind)                :: dold
        real(rkind)                :: dnew
        real(rkind)                :: mu
        real(rkind)                :: coeff
        real(rkind)                :: growth
        real(rkind)                :: decay
        real(rkind)                :: rk(3)

        real(rkind), parameter     :: pi = 3.14159265
        real(rkind), parameter     :: alpha = 3.0          ! Coefficient
        real(rkind), parameter     :: p  = 1.0
        real(rkind), parameter     :: q  = 0.5 
        real(rkind), parameter     :: Fc  = 2.0      ! Characteristic fractal dimension ! formel4 seite 59 khelifa and hill 
        real(rkind), parameter     :: dfc = 2000.E-6    ! Characteristic size of floc ! 2microMeter laut formel4 seite 59 khelifa and hill
        real(rkind), parameter     :: Fy  = 1.E-10   ! 10^-10 Yield strength of floc ( N )
        real(rkind), parameter     :: ka = 0.98            ! emp. coeff.
        real(rkind), parameter     :: kb = 3.3E-5      ! emp. coeff.
        real(rkind), parameter     :: ka2 = 14.6!0.98            ! emp. coeff.
        real(rkind), parameter     :: kb2 = 14.E3!3.3E-5      ! emp. coeff.
        real(rkind), parameter     :: dequi = 300E-6

        integer                    :: il

        g_dot_thr = 0.0000001
!AR: take care with dfv should be mu

        cnpar = 1

        do il=kbp(ip),nvrt

           dh = abs(znl(il+1,ip)-znl(il,ip))

           g2_dot = du + dv
           g_dot = sqrt(g2_dot)

           dold = dmf_mud(1,il,ip)

!            g_dot = 0.2     !DS: wasn hier los

           if (dmf_mud(1,il,ip) .lt. 0.) then
             write(errmsg,'(a20,I10,10F20.10)')'neg. diam in dmf_mud rk0', &
     & il, g_dot, beta, dm0, dfv(ip,il), dt, dold, dnew,nf(il), dold
             call parallel_abort(errmsg)
           endif

           if (dmf_mud(1,il,ip) .ne. dmf_mud(1,il,ip) ) then
             write(errmsg,'(a20,I10,10F15.10)')'NaN in dmf_mud rk0', il, nf(il), dold/dm0, g_dot, beta, dm0, dfv(ip,il), dt, dold, dnew
             call parallel_abort(errmsg)
           endif

           mu = 1.E-3  !*(prho(l,ip)+1000.)  !DS: was ist hier los
           beta = log(Fc/3.d0)/log(dfc/dm0)

           nf(il)= alpha * (dold/dm0)**beta
           coeff=(g_dot*dm0**beta)/(beta*log(dold/dm0)+1.d0)
           growth=0.3333*ka*mudc(1,il,ip)/rhosed*dm0**(nf(il)-3)*dold**(-nf(il)+4.-beta)
           decay=0.3333*kb*(mu*g_dot/fy)**0.5*dm0**(-1)*dold**(-beta+2)*(dold-dm0)
!           dddt1 = sign(min(0.1*dequi,abs(coeff*(growth-decay))),      &
!     &             coeff*(growth-decay)) 
           dddt1 = coeff*(growth-decay)
           dnew = dold + dddt1 * dt
           dold = dnew
           !dmf_mud(1,l,ip) = dnew

           if (dnew .ne. dnew ) then
             write(*,'(a20,I10,10F15.10)')'NaN in dmf_mud rk1', il, nf(il), dold/dm0, g_dot, beta, dm0, dfv(ip,il), dt,dold,dnew
             call parallel_abort(errmsg)
           endif

           if (dnew .lt. 0.) then
             write(errmsg,'(a20,I10,10F20.10)')'neg. diam in dmf_mud rk1', il, g_dot, beta, dm0, dfv(ip,il), dt, dold, dnew, nf(il), dold
             call parallel_abort(errmsg)
           endif

           nf(il)= alpha * (dold/dm0)**beta
           coeff=(g_dot*dm0**beta)/(beta*log(dold/dm0)+1.d0)
           growth=0.3333*ka*mudc(1,ip,il)/rhosed*dm0**(nf(il)-3)*dold**(-nf(il)+4.-beta)
           decay=0.3333*kb*(mu*g_dot/fy)**0.5*dm0**(-1)*dold**(-beta+2)*(dold-dm0)
!           dddt2 = sign(min(0.1*dequi,abs(coeff*(growth-decay))),      & 
!     &             coeff*(growth-decay))
           dddt2 = coeff*(growth-decay)
           dnew = 0.75*dmf_mud(1,ip,il)+0.25*dnew+0.25*dddt2*dt
           dold = dnew

           if (dnew .ne. dnew) then
             write(errmsg,'(a20,I10,10F15.10)')'NaN in dmf_mud rk2', il, nf(il), dold/dm0, g_dot, beta, dm0, dfv(ip,il), dt, dold, dnew
             call parallel_abort(errmsg)
           endif

           if (dnew .lt. 0.) then
             write(errmsg,'(a20,I10,10F20.10)')'neg. diam in dmf_mud rk2', il, g_dot, beta, dm0, dfv(ip,il), dt, dold, dnew, nf(il), dold
             call parallel_abort(errmsg)
           endif

           nf(il)= alpha * (dold/dm0)**beta
           coeff=(g_dot*dm0**beta)/(beta*log(dold/dm0)+1.d0)
           growth=0.3333*ka*mudc(1,ip,il)/rhosed*dm0**(nf(il)-3)*dold**(-nf(il)+4.-beta)
           decay=0.3333*kb*(mu*g_dot/fy)**0.5*dm0**(-1)*dold**(-beta+2)*(dold-dm0)
!           dddt3 = sign(min(0.1*dequi,abs(coeff*(growth-decay))),      &
!     &             coeff*(growth-decay))
           dddt3 = coeff*(growth-decay)
           dnew = 0.3333*dmf_mud(1,il,ip)+0.6666*dnew+0.6666*dddt3*dt
           dmf_mud(1,il,ip) = dnew
           if (dmf_mud(1,il,ip) .ne. dmf_mud(1,il,ip)) then
             write(errmsg,'(a20,I10,10F20.10)')'NaN in dmf_mud rk3', il, g_dot, beta, dm0, dfv(ip,il), dt, dold, dnew, nf(il), dold
             call parallel_abort(errmsg) 
           endif

           if (dmf_mud(1,il,ip) .lt. 0.) then
             write(errmsg,'(a20,I10,10F20.10)')'neg. diam in dmf_mud rk3', il, g_dot, beta, dm0, dfv(ip,il), dt, dold, dnew, nf(il), dold
             call parallel_abort(errmsg)
           endif


           if (ldebug .and. iplg(ip) == testnode) then
             write(1116,'(A10,I10,10F20.10)') 'get_floc', il,            & 
     &                  dddt2*1.E6,dm0*1.E6,dold*1.E6,dnew*1.E6,        &
     &                  growth*1.E6,decay*1.E6,(growth-decay)*1.E6,     &
     &                  beta, nf(il), g_dot
           endif

        end do

        nf(kbp(ip)) = nf(kbp(ip)+1)
        dmf_mud(1,kbp(ip),ip) = dmf_mud(1,kbp(ip)+1,ip)

      end subroutine get_floc_diam
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine RISK(Ustar,R)
        use elfe_glbl, only : rkind
        implicit none

        real(rkind), intent(in)  :: ustar
        real(rkind), intent(out) :: R

        IF ( Ustar >1.5 ) then
          R = 1.
        elseif ( Ustar <=0.7 ) then
          R = 0.
        ELSE
!         R = 1./(10.*Ustar**(-18)+1.)
          R = SQRT((Ustar-0.7)/0.8)
        endIF

      end subroutine RISK
!**********************************************************************
!*                                                                    *
!********************************************************************** 
      subroutine shields(ip,il,d,nu,rhost,usc)
        use elfe_glbl, only : rkind, grav, prho
        use flmud_pool !, only : rhow 
        implicit none

        integer, intent(in)      :: ip, il
        real(rkind), intent(in)  :: nu, d
        real(rkind), intent(out) :: rhost, usc

        real(rkind)              :: rgn , frstc, rhosusp, dst

        rhost=(prho(il,ip)-rhow)/rhow
        rgn = (Rhost*grav/nu**2)**(1./3.)
        dst = rgn*d
        if ( dst<6. ) then
          frstc = .109*dst**(-.5)
        elseif ( dst<10. ) then
          frstc = .14*dst**(-.64)
        elseif ( dst<20. ) then
          frstc = .04*dst**(-.1)
        elseif ( dst<150. ) then
          frstc = .013*dst**(.29)
        elseif ( dst>=150. ) then
          frstc = .055
        endIF
        usc = SQRT(frstc*(rhost*grav*d))

      end subroutine shields
!**********************************************************************
!*                                                                    *
!**********************************************************************
       subroutine get_cgel(ip,nf,cgel,phigel)
        use elfe_glbl, only : rkind, grav, kbp, nvrt,iplg
        use flmud_pool
        implicit none

        real(rkind), intent(in)  :: nf(nvrt)    !Fractal dimensionr
        real(rkind), intent(out) :: cgel(nvrt)  ! gel, concentration
        real(rkind), intent(out) :: phigel(nvrt)! gel, concentration

        integer, intent(in) :: ip
!        integer, intent(in) :: testnode
!        logical, intent(in) :: ldebug

        integer :: il

        do il = kbp(ip), nvrt
          !cgel(il) = rhosed * (dm0/dmf_mud(1,il,ip))**(3.-nf(il))
          cgel(il) = rhomud(1,il,ip)* (dm0/dmf_mud(1,il,ip))**(3.-nf(il))
          phigel(il) = 1.
          if (ldebug .and. iplg(ip) == testnode) then
          write(1113,'(A10,I10,8F15.8)') 'get_cgel',il,  &
     &            rhosed,dmf_mud(1,il,ip),dm0,nf(il),      &
     &           (dmf_mud(1,il,ip)/dm0)**(3.-nf(il)),      &
     &            cgel(il),phigel(il)
           endif
        end do

       end subroutine get_cgel
!**********************************************************************
!*                                                                    *
!**********************************************************************
        subroutine set_yield(ip,rough_p,tstress)
        use elfe_glbl, only : rkind, nvrt, kbp, dfv, dfh, prho,iplg
        use flmud_pool !, only : dfv_yield, dfh_yield, rough_mud
        implicit none

        integer, intent(in) :: ip
!        integer, intent(in) :: testnode
!        logical, intent(in) :: ldebug

        real(rkind), intent(in)   :: rough_p !(nvrt)
        real(rkind), intent(in)   :: tstress(nvrt)

        real(rkind)         :: smooth, rhop, tau0, fak

        integer             :: il 

        smooth = 0.000001

        do il=kbp(ip),nvrt
          rhop = max(0.d0,dble(prho(il,ip)))/1000.
          call set_yieldstress(rhop,tau0)
          if (tau0 .gt. 0. .and. rhop .gt. 0.) then
            fak = 0.5d0+0.5d0*tanh((tstress(ip)-tau0)/smooth)
            dfv_yield(ip,il) = fak * dfv(ip,il) + (1.-fak) * 99.
            dfh_yield(ip,il) = fak * dfh(ip,il) + (1.-fak) *  0.
            if (il==kbp(ip)) rough_mud(ip) = fak * rough_p + (1.-fak) * 99.
          else
            dfv_yield(ip,il) = dfv(ip,il)
            dfh_yield(ip,il) = dfh(ip,il)
            if (il==kbp(ip)) rough_mud(ip) = rough_p !(ip)
          endif
          if (ldebug .and. iplg(ip) == testnode) then
            write(1117,'(A10,I5,17F14.8)') 'set yield',  &
     &      il, rhop,tstress(ip)-tau0,tstress(ip),tau0, &
     &      dfv_yield(ip,il),dfh_yield(ip,il), &
     &      dfv(ip,il),dfh(ip,il),rough_p,rough_mud(ip),fak, &
     &      tanh((tstress(ip)-tau0)/smooth)
          endif
        end do

      end subroutine
!**********************************************************************
!*                                                                    *
!**********************************************************************
        subroutine write_point_output(ip,SS1d,NN1d,tke1d,eps1d,L1d,num1d,nuh1d)
        use elfe_glbl, only : rkind, nvrt, kbp, dfv, dfh, prho,iplg,uu2,vv2,znl,wsink,vts
        use flmud_pool 
        implicit none

        integer, intent(in) :: ip
        real(rkind), intent(in), dimension(0:nvrt) :: SS1d,NN1d,tke1d,eps1d,L1d,num1d,nuh1d
        real(rkind) :: rstress,rich,rif
        integer             :: i,i2
        integer, save       :: iwrite
        data iwrite/0/

        if (iplg(ip)==testnode.and.ldebug) then
!          iwrite = iwrite + 1
!          if (mod(iwrite,icycle) == 0) then
            do i=kbp(ip),nvrt
              i2=i-kbp(ip) !from GOTM
              rstress=sqrt(SS1d(i2))*dfv(ip,i) !Reynolds stress
              if(SS1d(i2)==0) then
                rich=0
              else
                rich=NN1d(i2)/sqrt(SS1d(i2)) 
              endif
              rif=0 !flux Rich.
              write(1020,'(I10,400(1x,e15.7))') i,znl(i,ip),tke1d(i2),eps1d(i2),num1d(i2), &
     &                     nuh1d(i2),prho(i,ip),-wsink(1,i,ip)*1000,vts(i,ip),rstress, &
     &                     mudc(1,i,ip),-uu2(i,ip),dmf_mud(1,i,ip)*1000,rich,rif
            end do
!          endif !mod
        endif !ip
      end subroutine
!**********************************************************************
!*                                                                    *
!**********************************************************************

