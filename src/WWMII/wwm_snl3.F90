#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine triadswan (ip, hs, smespc, acloc, imatra, imatda, ssnl3)
!
      use datapool
      implicit none

      integer, intent(in)        :: ip
      real(rkind), intent(in)    :: hs, smespc
      real(rkind), intent(out)   :: ssnl3(msc,mdc)
      real(rkind), intent(in)    :: acloc(msc,mdc)
      real(rkind), intent(inout) :: imatra(msc,mdc), imatda(msc,mdc)
!
      integer i1, i2, id, is, ism, ism1, ismax, isp, isp1,  ij1, ij2, ires

      real(rkind)    aux1, aux2, biph, c0, cm, dep_2, dep_3, e0
      real(rkind)    em,ft, rint, sigpi, sinbph, stri, wism, wism1 , fac1
      real(rkind)    wisp, wisp1,w0, wm, wn0, wnm,  xisln, ursell
      
      real(rkind), allocatable :: e(:), sa(:,:)

      ptriad(1)  = 0.05_rkind    ! swan settings 40.51
      ptriad(2)  = 2.5_rkind     ! frequency range for triad = trira * smebrk = ismax
      ptriad(3)  = 10._rkind     ! not used
      ptriad(4)  = 0.2_rkind     ! parameterized biphase
      ptriad(5)  = 0.01_rkind    ! ursell lower limit

      if (trico .gt. 0.)   ptriad(1) = trico
      if (trira  .gt. 0.)  ptriad(2) = trira
      if (triurs .gt. 0.)  ptriad(5) = triurs

      call ursell_number(hs,smespc,dep(ip),ursell)

      if (ursell .lt. ptriad(5)) return

      !write(*,'(5f15.6)') dep(ip), smespc, ursell, ptriad(5)

      dep_2 = dep(ip)**2
      dep_3 = dep(ip)**3
      i2     = int (float(msc) / two)
      i1     = i2 - 1
      xis    = spsig(i2) / spsig(i1)
      xisln  = log( xis )
      isp    = int( log(two) / xisln )
      isp1   = isp + 1
      wisp   = (two - xis**isp) / (xis**isp1 - xis**isp)
      wisp1  = one - wisp
      ism    = int( log(0.5_rkind) / xisln )
      ism1   = ism - 1
      wism   = (xis**ism -0.5_rkind) / (xis**ism - xis**ism1)
      wism1  = 1. - wism

      allocate (e (1:msc))
      allocate (sa(1:msc+isp1,1:mdc))
      e  = zero 
      sa = zero 

      ismax = 1
      do is = 1, msc
       if ( spsig(is) .lt. ( ptriad(2) * smespc) ) then
          ismax = is
        endif
      enddo
!
!      ismax = max( 1, min( msc, max ( ismax , isp1 ) ) ) ! added fix the bug described below ...
      ismax = max ( ismax , isp1 ) 
!
        biph   = (0.5_rkind*pi)*(tanh(ptriad(4)/ursell)-1.)
        sinbph = abs( sin(biph) )
        do id = 1, mdc
           do is = 1, msc
              e(is) = acloc(is,id) * pi2 * spsig(is)
           end do
           do is = 1, ismax 
              e0  = e(is)
              w0  = spsig(is)
              wn0 = wk(ip,is)
              c0  = w0 / wn0
              if ( is.gt.-ism1 ) then
                 em  = wism * e(is+ism1)      + wism1 * e(is+ism)
                 wm  = wism * spsig(is+ism1)  + wism1 * spsig(is+ism)
                 wnm = wism * wk(ip,is+ism1)  + wism1 * wk(ip,is+ism)
                 cm  = wm / wnm
              else
                 em  = zero 
                 wm  = zero 
                 wnm = zero 
                 cm  = zero 
              end if
              aux1 = wnm**2 * ( g9 * dep(ip) + two*cm**2 )
              aux2 = wn0 * dep(ip) * ( g9 * dep(ip) + (two/15._rkind) * g9 * dep_3 * wn0**2 -(two/five) * w0**2 * dep_2 )
              rint = aux1 / aux2
              ft = ptriad(1) * c0 * cg(ip,is) * rint**2 * sinbph
              sa(is,id) = max(zero, ft * ( em * em - two * em * e0 ))
           end do
        end do

        do is = 1, msc
           sigpi = spsig(is) * pi2
           do id = 1, mdc
             if (acloc(is,id) .gt. verysmall) then
               stri = sa(is,id) - 2.*(wisp  * sa(is+isp1,id) + wisp1 * sa(is+isp,id))
               if (abs(stri) .gt. 0.) then
                 if (icomp .ge. 2) then
                   if (stri .gt. 0.) then
                     imatra(is,id) = imatra(is,id) + stri / sigpi
                   else
                     imatda(is,id) = imatda(is,id) - stri / (acloc(is,id)*sigpi)
                   end if
                 else
                   imatra(is,id) = imatra(is,id) + stri / sigpi
                   imatda(is,id) = imatda(is,id) + stri / (acloc(is,id)*sigpi)
                 end if
                 ssnl3(is,id) = stri/sigpi
               end if
            end if
          end do
        end do

      deallocate(e,sa)

      end subroutine 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine triad_dingemans (ip, acloc, imatra, imatda, ssnl3)
!
        use datapool
        implicit none

        integer, intent(in) :: ip
        real(rkind), intent(in)    :: acloc(msc,mdc)
        real(rkind), intent(inout) :: imatra(msc,mdc), imatda(msc,mdc)
        real(rkind), intent(out)   :: ssnl3(msc,mdc)
        integer             :: is, is2, id
        real(rkind)         :: ecloc(msc,mdc), e2(msc,mdc), d20
        real(rkind)         :: df, domega, omega, omega1, fac, z1a, z1b
        integer             :: j, j1, j2, jmin, j2abs

        do is = 1, msc
          do id = 1, mdc 
            ecloc(is,id) = acloc(is,id) * spsig(is) * ddir
          end do 
        end do

        ssnl3 = 0.

        do id = 1, mdc
          do is = 1,msc-1
            df = ((spsig(is+1) - spsig(is)))/pi2
            domega = pi2 * df
            omega = spsig(is) * pi2 
            jmin = nint(0.5*MyREAL(is))
            if (2*jmin .eq. is) then
               fac = 0.5 
            else
               fac = 1. 
            endif
            z1a = dep(ip) * (jmin*domega)**2 / g9
            z1b = dep(ip) * ((jmin-is-1)*domega)**2 / g9
            e2(is,id) = 0. 
            do is2 = jmin, msc
               j1    = is2 
               j2    = is2 - is 
               j2abs = iabs (j2)
!            zero frequencies are skipped
               if (j2 .eq. 0)  cycle
               omega1 = is2 * domega
               e2(is,id) = e2(is,id) + fac * d20(omega,omega1,dep(ip),z1a,z1b)* ecloc(j1,id)*ecloc(j2abs,id)
               fac    = 1. 
            end do
            e2(is,id) = e2(is,id) * df
          end do
        end do
      end subroutine
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function d20(omega, omega1, h, z1a, z1b)
      use datapool, only : g9, rkind
      implicit none
!
!----------------------------------------------------------------------*
!     compute the second order spectral response function.             *
!--------------------------------- delft hydraulics -- gkl -- 880419 --*
!
      real(rkind), intent(in) :: omega, omega1, h
      real(rkind), intent(inout) :: z1a, z1b
      real(rkind) :: k, k1, k2
      real(rkind) :: omega2
      real(rkind) :: aome2
      real(rkind) :: cthkh, cthk1h, cthk2h 
      real(rkind) :: d2, d20
!
      omega2 = omega - omega1
      aome2  = abs (omega2)
      z1a = abs(z1a)
      z1b = abs(z1b)
      call dispu2 (omega1, h, z1a, k1)
      call dispu2 (aome2,  h, z1b, k2)
      if (omega2 .lt. 0.) k2 = - k2
      k = k1 + k2
      cthk1h = 1. / tanh (k1*h)
      cthk2h = 1. / tanh (k2*h)
      cthkh  = 1. / tanh (k*h)
      d2 = 0.5 *  (omega1**2 + omega2**2 + omega1*omega2 -              &
     &             omega1 * omega2 * cthk1h * cthk2h - omega *          & 
     &             omega1 * cthkh  * cthk1h - omega  * omega2 *         &
     &             cthkh  * cthk2h)
      d20 = d2 / (g9 * (1. - omega**2 * cthkh / (g9*k)))
      d20 = 2. * d20**2
!
      return
      end function 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine dispu2(omega,d,z1,k)
      use datapool, only : rkind
      implicit none

      real(rkind), intent(in)    :: omega, d
      real(rkind), intent(inout) :: z1
      real(rkind), intent(out)   :: k
      real(rkind), parameter     :: eps = 0.0001
      real(rkind), parameter     :: g = 9.81

      real(rkind) :: z0, z2, fak1, fak2, sig 

      z0 = d*omega*omega/g
   10    sig = tanh(z1)
         fak1 = z1*sig
         fak2 = z1 + sig*(1.-fak1)
         z2 = z1 + (z0-fak1)/fak2
         if (abs((z2-z1)/z2).gt.eps) goto 40
         goto 60
   40    z1 = z2
         goto 10
   60 k = z2/d
      z1 = z2
      return
      end subroutine
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function delta(ip, is, is1, is2) result(res)
      use datapool, only : rkind, wk
      implicit none

      integer, intent(in)        :: ip, is, is1, is2
      real(rkind)                :: res !return

        res = wk(ip,is) - wk(ip,is1) - wk(ip,is2)
        
      end function delta
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function ddelta_dx(ip, is, is1, is2, id) result(res)
      use datapool, only : rkind, sinth, costh, dwkdx, dwkdy
      implicit none

      integer, intent(in)        :: ip, is, is1, is2, id
      real(rkind)                :: res !return

      real(rkind)                :: ka_1abl, kb_1abl, kc_1abl

        ka_1abl = costh(id) * dwkdx(ip,is) + sinth(id) * dwkdy(ip,is)
        kb_1abl = costh(id) * dwkdx(ip,is1) + sinth(id) * dwkdy(ip,is1)
        kc_1abl = costh(id) * dwkdx(ip,is2) + sinth(id) * dwkdy(ip,is2)
      
        res = ka_1abl - kb_1abl - kc_1abl
        
      end function ddelta_dx
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function w(ip, is, is1, is2, n1, emf) result(res)
      use datapool, only : rkind, wk, cg, spsig, g9, one
      implicit none
     
      integer, intent(in)        :: ip, is, is1, is2, n1, emf
      real(rkind)                :: res !return

!!! mapping var.
      real(rkind)                :: omegaa  
      real(rkind)                :: omegab 
      real(rkind)                :: omegac 
      real(rkind)                :: cga  
      real(rkind)                :: cgb  
      real(rkind)                :: cgc  
      real(rkind)                :: ka  
      real(rkind)                :: kb  
      real(rkind)                :: kc 
!!! function declaration
      real(rkind)                :: tau
      
!       n1=0 means +
!       n1=1 means -
!       em=0 without eldeberky & madsen
!       em=1 with eldeberky & madsen

        omegaa = spsig(is)
        omegab = spsig(is)
        omegac = spsig(is)

        cga    = cg(ip,is)
        cgc    = cg(ip,is1)
        cgb    = cg(ip,is2)

        ka      = wk(ip,is)
        kb      = wk(ip,is1)
        kc      = wk(ip,is2)

!AR: reduce the +- stuff ...
!
        res =  (g9 * omegaa) / (8 * sqrt(cga * cgb *cgc))  *  &        
     &         ( (-1)**n1 * (2 - emf * tau(ip, is, is1, is2, n1)) * kb * kc +  & 
     &         ( 1 - emf * tau(ip, is, is1, is2, n1)) * ((omegab * omegac)**2 / g9**2) + & 
     &         ( kb**2 * omegac / omegaa) + ( (-1)**n1 * kc**2 * omegab / omegaa) + ( (-1)**(n1+1) * & 
     &         ( 1- emf * tau(ip, is, is1, is2,n1)) * ((omegaa**2 * omegab * omegac) / g9**2) ) )
        
      end function w
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function tau(ip, is, is1, is2, n1) result(res)
      use datapool, only : rkind, wk, cg, spsig
      implicit none

      integer, intent(in)        :: ip, is, is1, is2, n1
      real(rkind)                :: res !return

      real(rkind)                :: cg_a 
      real(rkind)                :: ka  
      real(rkind)                :: kb  
      real(rkind)                :: kc  
      real(rkind)                :: omega_a  

        omega_a = spsig(is)

        cg_a    = cg(ip,is)

        ka      = wk(ip,is)
        kb      = wk(ip,is1)
        kc      = wk(ip,is2)
      
        res = 2 * cg_a * ( ka + (-1)**(n1+1) * kb - kc ) / omega_a
            
      end function tau
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function dwdx(ip, is, is1, is2, id, n1, switch) result(res)
      use datapool, only : rkind, wk, cg, spsig, g9, one, costh, sinth, dcgdx, dcgdy, dwkdx, dwkdy, one, two, three
      implicit none
     
      integer, intent(in)        :: ip, is, is1, is2, id, n1
      real(rkind)                :: res !return

      integer                    :: switch

      real(rkind)                :: omegaa 
      real(rkind)                :: omegab 
      real(rkind)                :: omegac 
      real(rkind)                :: cga 
      real(rkind)                :: cgb 
      real(rkind)                :: cgc 
      real(rkind)                :: cga_
      real(rkind)                :: cgb_
      real(rkind)                :: cgc_
      real(rkind)                :: ka 
      real(rkind)                :: kb 
      real(rkind)                :: kc 
      real(rkind)                :: ka_ 
      real(rkind)                :: kb_ 
      real(rkind)                :: kc_ 

        omegaa = spsig(is)
        omegab = spsig(is)
        omegac = spsig(is)

        cga    = cg(ip,is)
        cgc    = cg(ip,is1)
        cgb    = cg(ip,is2)

        ka      = wk(ip,is)
        kb      = wk(ip,is1)
        kc      = wk(ip,is2)

        cga_ = costh(id) * dcgdx(ip,is) + sinth(id) * dcgdy(ip,is) 
        cgb_ = costh(id) * dcgdx(ip,is1) + sinth(id) * dcgdy(ip,is1)
        cgc_ = costh(id) * dcgdx(ip,is2) + sinth(id) * dcgdy(ip,is2)

        ka_ = costh(id) * dwkdx(ip,is) + sinth(id) * dwkdy(ip,is)
        kb_ = costh(id) * dwkdx(ip,is1) + sinth(id) * dwkdy(ip,is1)
        kc_ = costh(id) * dwkdx(ip,is2) + sinth(id) * dwkdy(ip,is2)

      
        res = (g9*((4*Cga*Cgb*Cgc*((-1)**n1*((kb_*kc + kb*kc_)*omegaa +kc*kc_*omegab) +kb*kb_*omegac))/&
               omegaa -(Cga*Cgb_*Cgc +Cgb*(Cga_*Cgc + Cga*Cgc_))*(2*(-1)**n1*kb*kc -((-1)**n1*omegaa**2*omegab*omegac)/g9**2+&
              (omegab**2*omegac**2)/g9**2 +((-1)**n1*kc**2*omegab +kb**2*omegac)/omegaa)))/(16.*(Cga*Cgb*Cgc)**1.5)

      end function dwdx        
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function k(ip, is, is1, is2, id, is3, is4, is5, n1, emf) result(res)
      use datapool, only : rkind, g9, one
      implicit none

      integer, intent(in)        :: ip, is, is1, is2, is3, is4, is5, id, n1, emf
      real(rkind)                :: res ! return
      real(rkind)                :: delta, dwdx, w, ddelta_dx
      
       res = one / ((delta(ip, is3, is4, is5))**2) * dwdx(ip,is,is1,is2,id,n1,emf) - (w(ip,is,is1,is2,n1,emf) / ((delta(ip, is, is1, is2))**3)) * ddelta_dx(ip, is3, is4, is5, id)

      end function k
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function j(ip, is, is1, is2, id) result(res)
      use datapool, only : rkind, one
      implicit none
     
      integer, intent(in)        :: ip, is, is1, is2, id 
      real(rkind)                :: res ! return 
      !!! function declaration 
      real(rkind)                :: delta, ddelta_dx
           
        res = - ( one / (delta(ip, is, is1, is2)**3)  ) * ddelta_dx(ip, is, is1, is2, id)
        
      end function j
!**********************************************************************
!*                                                                    *
!**********************************************************************
     subroutine snl3ta(ip,snl3,dsnl3)
     use datapool, only : rkind, msc, mdc, ac2, ZERO, spsig, cg, frintf, ddir
     implicit none

     real(rkind), intent(out) :: snl3(msc,mdc), dsnl3(msc,mdc)
     integer, intent(in)      :: ip
     integer :: is, is1, is2, id, em, kron_delta
     real(rkind) :: SUPER, SUB, f, f1, f2, tt, j, SUPERD, SUBD, k, w, JAC



     do id = 1, mdc
       snl3(:,id) = zero
       dsnl3(:,id) = zero
       SUPER = ZERO; SUPERD = ZERO
       SUB   = ZERO; SUBD = ZERO
       do is = 1, msc
         f = (ac2(ip,is,id) * spsig(is)**2. * frintf * ddir)* cg(ip,is)
         do is1 = 1, msc
           f1 = (ac2(ip,is1,id) * spsig(is1)**2. * frintf * ddir) * cg(ip,is1)
           do is2 = 1, msc
             f2 = (ac2(ip,is2,id) * spsig(is2)**2. * frintf * ddir) * cg(ip,is2)
             SUPER = SUPER + ( k(ip, is, is1, is2, id, is, is1, is2, 0, 0) * f1 * f2 + k(ip, is1, is2, is, id, is, is1, is2, 1, 0) * f1 * f &
     &                       + k(ip, is2, is1, is, id, is, is1, is2, 1, 0) * f2 * f ) * w(ip, is, is1, is2, 0, 0) * kron_delta(is,is1+is2)
             SUPERD = SUPERD +  ( k(ip, is1, is2, is, id, is, is1, is2, 1, 0) * f1 + &
     &                            k(ip, is2, is1, is, id, is, is1, is2, 1, 0) * f2 ) * w(ip, is, is1, is2, 0, 0) * kron_delta(is,is1+is2) 
           enddo
         enddo
         SUPER = 4 * SUPER
         do is1 = 1, msc
           f1 = ac2(ip,is1,id)
           do is2 = 1, msc
             f1 = ac2(ip,is2,id)
             SUB = SUB + ( k(ip, is, is1, is2, id, is2, is, is1, 1, 0) * f1 * f2 + k(ip, is1, is, is2, id, is2, is, is1, 1, 0) * f2 * f + &
      &                    k(ip, is2, is1, is, id, is2, is, is1, 0, 0) * f1 * f ) * w(ip, is, is1, is2, 1, 0) * kron_delta(is2,is+is1)
             SUBD = SUBD + ( k(ip, is, is1, is2, id, is2, is, is1, 1, 0) * f1 * f2 + k(ip, is1, is, is2, id, is2, is, is1, 1, 0) * f2 * f + &
      &                      k(ip, is2, is1, is, id, is2, is, is1, 0, 0) * f1 * f ) * w(ip, is, is1, is2, 1, 0) * kron_delta(is2,is+is1)
           enddo
         enddo
         SUB = 8 * SUB
         JAC = 1./(SPSIG(IS)*DDIR*FRINTF)
         snl3(is,id) = (SUB + SUPER) * JAC
         dsnl3(is,id) = (SUB + SUPER) 
       enddo ! all freq. 
     enddo ! all directions 

     end subroutine snl3ta
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function kron_delta(i, j) result(res)
      implicit none

      integer, intent(in)        :: i,j
      integer                    :: res ! return 

        res = int((float(i+j)-abs(i-j)))/(float((i+j)+abs(i-j))) 

      end function kron_delta 
!**********************************************************************
!*                                                                    *
!**********************************************************************
