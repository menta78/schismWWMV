      SUBROUTINE sed_bedload_mpm()

! this whole section is currently unused
! commented for better overview
!... Compute critical stress for horizontal bed
!... For Meyer-Peter Muller formulation tauc0= 0.047.
!... Compute bottom stress and tauc. Effect of bottom slope (Carmo,1995).
!            tauc0 =tau_ce(ised)
!            alphas=datan(-(dzdx*angleu+dzdy*anglev))

!... Magnitude of bed load at rho points. Meyer-Peter Muller formulation.
!... bedld has dimensions  (m2 s-1)
!!!#if defined CARMO
!!!          if(slope_formulation==sf_carmo)
!!!            call stress(angleu,anglev,dzdx,dzdy,alphas,sed_angle,tauc0,tauc,i)
!!!            Eq. 46 q_{b,q} = 8*(tau_k
!!!            bedld=8.0_r8*(MAX((tau_w(i)-tauc)/smgd,0.0_r8)**1.5_r8)*smgdr
!!!#elif defined SOULSBY
!!!          elseif(slope_formulation==sf_soul)
!!!            call stress_soulsby(angleu,anglev,dzdx,dzdy,alphas,sed_angle,tauc0,tauc,i)
!!!            bedld=8.0_r8*(MAX((tau_w(i)-tauc)/smgd,0.0_r8)**1.5_r8)*smgdr
!!!#elif defined DELFT
!!!          elseif(slope_formulation==sf_del)
!!!            bedld=8.0_r8*(MAX((tau_w(i)*osmgd-0.047_r8),0.0_r8)**1.5_r8)*smgdr!!!          endif

!!!#elif defined DAMGAARD
!!!            if (abs(alphas)>datan(sed_angle))alphas=datan(sed_angle)*SIGN(1.0_r8,alphas)
!!!            tauc=tauc0*(sin(datan(sed_angle)+(alphas))/sin(datan(sed_angle)))
!!!            bedld=8.0_r8*(MAX((tau_w(i)-tauc)/smgd,0.0_r8)**1.5_r8)*smgdr

!!!            if (alphas>0)then
!!!                 cff=1
!!!             else if (alphas<=0)then
!!!                 cff=1+0.8*((tauc0/tau_w(i))**0.2)*(1-(tauc/tauc0))**(1.5+(tau_w(i)/tauc0))
!!!             endif
!!!             bedld=bedld*cff
!!!             if (time>1500) then
!!!!ZYL
!!!               if (isnan(cff)==.true.) call parallel_abort('cff is NaN "i" 1')
!!!               if (isnan(bedld)==.true.) call parallel_abort('bedld is NaN "i"')
!!!!"
!!!             endif
!!!#endif  ! bedld

!!!!-------------------------------------------------------------------
!!!!... Partition bedld into x  and y directions, at the center of each 
!!!!... element (rho points), and integrate in time.
!!!!... FX_r and FY_r have dimensions of m2
!!!!-------------------------------------------------------------------

!!!#if defined CARMO || defined SOULSBY
!!!          if(slope_formulation==3) ! || soulsby
!!!            FX_r(i)=bedld*dcos(alphas)*angleu*dt
!!!            FY_r(i)=bedld*dcos(alphas)*anglev*dt
!!!#elif defined DAMGAARD || defined DELFT
!!!          else
!!!            FX_r(i)=bedld*angleu*dt
!!!            FY_r(i)=bedld*anglev*dt
!!!          endif
!!!#endif !CARMO||SOULSBY

!!!#ifdef DELFT
!!!          if(slope_formulation== 
!!!!... Bed_slope effects
!!!!... longitudinal bed slope
!!!!... limit slope to 0.9*(sed_angle)

!!!            cff=(dzdx*angleu+dzdy*anglev)
!!!            cff1=min(abs(cff),0.9*sed_angle)*sign(1.0_r8,cff)

!!!!... add contribution of longitudinal bed slope to bed load

!!!            cff2=datan(cff1)
!!!            a_slopex=1+1*((sed_angle/(cos(cff2)*(sed_angle-cff1)))-1)

!!!            FX_r(i)=FX_r(i)*a_slopex
!!!            FY_r(i)=FY_r(i)*a_slopex

!!!!... Transverse bed slope

!!!            cff=(-(dzdx*anglev)+dzdy*angleu)
!!!            a_slopey=1.5*sqrt(tauc0/(abs(bustr(i))+abs(bvstr(i))))*cff
!!!            FX_r(i)=FX_r(i)-(FY_r(i)*a_slopey)
!!!            FY_r(i)=FY_r(i)+(FX_r(i)*a_slopey)
!!!#endif !DELFT
      END SUBROUTINE sed_bedload_mpm
