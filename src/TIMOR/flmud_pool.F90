!**********************************************************************
!*                                                                    *
!**********************************************************************
      MODULE flmud_pool
        use elfe_glbl, only : rkind,mudc=>tr_nd,rhosed,laddmud_d,laddmud_v,vts,wsink,rhomud

        save

        integer                          :: ithick !# of bottom levels for init. & dump mud
        logical                          :: lhindws !sink vel. with hindrance (otherwise Stokes law only)
        logical                          :: ldumpmud
        logical                          :: lfloc !call floc routine
        logical                          :: ldebug
        integer                          :: testnode
!        integer                          :: icycle
!        logical                          :: laddmud_d !switch on/off density effects
!        logical                          :: laddmud_v !switch on/off viscosity effects

        real(rkind)                      :: dm0 !primary floc diameter
        real(rkind)                      :: dm1 !initial particle diameter
        real(rkind)                      :: mudref !ref. mud conc [kg/m3]

!        real(rkind), parameter     :: rhosed = 2650.   ! kg/m3
        real(rkind), parameter     :: rhow   = 1000.
!        real(rkind), allocatable   :: rhomud(:,:,:)     ! rhomud(nlvdim,nkndim)  !Mud floc particle density [kg/m3]
        real(rkind), allocatable   :: dmf_mud(:,:,:)    ! dmf_mud(ntracers,nvrt,npa) !floc diameter [m]
        real(rkind), allocatable   :: nfg(:,:)  !Fractal dimension
!       mudc(ntracers,nvrt,npa)   !Mud concentration [kg/m^3]
!        real(rkind), allocatable   :: wsink(:,:,:)      ! wsink(nlvdim,nkndim)   !Sink velocity [m/s]
!        real(rkind), allocatable   :: vts(:,:)        ! vts(nlvdim,nkndim)     !rheological viscosity [m^2/s]
        real(rkind), allocatable   :: btmshear(:)     ! bottom shear stress    !bottom shear stress [Pa]
        real(rkind), allocatable   :: rough_mud(:)    ! roughness due to mud   !roughness length [m]
        real(rkind), allocatable   :: dfv_yield(:,:)  ! visc. in yield-regime  
        real(rkind), allocatable   :: dfh_yield(:,:)  ! diff. in yield-regime  

       END MODULE flmud_pool
!**********************************************************************
!*                                                                    *
!**********************************************************************
