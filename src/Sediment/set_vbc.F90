      SUBROUTINE set_vbc()
!--------------------------------------------------------------------!
! This routine sets vertical boundary conditons for tracers.         !
!                                                                    !
! This subroutine is adapted from a ROMS routine                     !
! Copyright (c) 2002-2007 The ROMS/TOMS Group                        !
!   Licensed under a MIT/X style license                             !
!   See License_ROMS.txt                                             !
!                                                                    !
! Author: Ligia Pinto                                                !
! Date: xx/08/2007                                                   !
!                                                                    !
! History: 2012/12 - F.Ganthy : form homogenisation of sediments     !
!          routines                                                  !
!                                                                    !
!--------------------------------------------------------------------!

      USE sed_mod
      USE elfe_glbl, ONLY: rkind,nvrt,nea,dfv,idry_e,kbe,nm,uu2,vv2,ze
      USE elfe_msgp, ONLY: myrank,parallel_abort

      IMPLICIT NONE
      SAVE

!- Local variables --------------------------------------------------!

      INTEGER     :: i,j,n1,n2,n3
      REAL(rkind) :: cff1, cff2, cff3, cff4, cff5,tmp
      REAL(rkind) :: wrk

!- Start Statement --------------------------------------------------!


!---------------------------------------------------------------------
! - Set kinematic bottom momentum flux (m2/s2).
!---------------------------------------------------------------------

      IF (drag_formulation == 1) THEN
      ! - Set logarithmic bottom stress.

        DO i=1,nea
          IF(idry_e(i)==1) CYCLE

          tmp = ze(kbe(i)+1,i)-ze(kbe(i),i)                

          IF(tmp>Zob(i)) THEN
            IF(tmp<=0.OR.Zob(i)<=0) THEN
              CALL parallel_abort('SEDIMENT: Cd failed')
            ENDIF
            cff1 = 1.0d0/DLOG(tmp/Zob(i))
            cff2 = vonKar*vonKar*cff1*cff1
            wrk  = MIN(Cdb_max,MAX(Cdb_min,cff2))
          ENDIF

          n1=nm(i,1)
          n2=nm(i,2)
          n3=nm(i,3)

          IF((ze(kbe(i)+1,i)-ze(kbe(i),i))>Zob(i)) THEN
            cff3 = (vv2(kbe(i)+1,n1)+                                &
            &       vv2(kbe(i)+1,n2)+                                &
            &       vv2(kbe(i)+1,n3))/3.0d0

            cff4 = (uu2(kbe(i)+1,n1)+                                &
            &       uu2(kbe(i)+1,n2)+                                &
            &       uu2(kbe(i)+1,n3))/3.0d0

            cff5 = SQRT(cff4*cff4+cff3*cff3)

            ! Bottom stress
            bustr(i) = wrk*cff4*cff5
            bvstr(i) = wrk*cff3*cff5

          ELSE

            cff3 = (vv2(kbe(i)+2,n1)+                                &
            &       vv2(kbe(i)+2,n2)+                                &
            &       vv2(kbe(i)+2,n3))/3.0d0-                         &
            &      (vv2(kbe(i)+1,n1)+                                &
            &       vv2(kbe(i)+1,n2)+                                &
            &       vv2(kbe(i)+1,n3))/3.0d0

            cff4 = (uu2(kbe(i)+2,n1)+                                &
            &       uu2(kbe(i)+2,n2)+                                &
            &       uu2(kbe(i)+2,n3))/3.0d0-                         &
            &      (uu2(kbe(i)+1,n1)+                                &
            &       uu2(kbe(i)+1,n2)+                                &
            &       uu2(kbe(i)+1,n3))/3.0d0

            cff5 = (dfv(n1,kbe(i)+1)+                                &
            &       dfv(n2,kbe(i)+1)+                                &
            &       dfv(n3,kbe(i)+1)+                                &
            &       dfv(n1,kbe(i)+2)+                                &
            &       dfv(n2,kbe(i)+2)+                                &
            &       dfv(n3,kbe(i)+2))/6.0d0

            tmp = ze(kbe(i)+2,i)-ze(kbe(i)+1,i)

            IF(tmp<=0) CALL parallel_abort('SEDIMENT: div. by 0')
            ! Bottom stress
            bustr(i) = cff5*cff4/tmp
            bvstr(i) = cff5*cff3/tmp

          ENDIF ! End test on Zob>ze

        END DO !End loop nea

!      ELSEIF (drag_formulation == 2) THEN
!      ! Set quadratic bottom stress.
!
!        DO i=1,nea
!          if(idry_e(i)==1) cycle
!          n1=nm(i,1)
!          n2=nm(i,2)
!          n3=nm(i,3)
!          cff1=(vv2(kbe(i)+1,n1)+vv2(kbe(i)+1,n2)+vv2(kbe(i)+1,n3))/3
!          cff2=(uu2(kbe(i)+1,n1)+uu2(kbe(i)+1,n2)+uu2(kbe(i)+1,n3))/3

!          cff3=SQRT(cff2*cff2+cff1*cff1)
!!Error: rdrg2 not defined
!!LLP user defined, if defined valor must be read from ...
!         bustr(i)=rdrg2*cff2*cff3
!         bvstr(i)=rdrg2*cff1*cff3
!       END DO

!       else if (drag_formulation == 3) then
!!#elif defined UV_LDRAG
!!
!!  Set linear bottom stress.
!!
!       DO i=1,nea
!         if(idry_e(i)==1) cycle
!         n1=nm(i,1)
!         n2=nm(i,2)
!         n3=nm(i,3)
!         cff1=(vv2(kbe(i)+1,n1)+vv2(kbe(i)+1,n2)+vv2(kbe(i)+1,n3)/3
!         cff2=(uu2(kbe(i)+1,n1)+uu2(kbe(i)+1,n2)+uu2(kbe(i)+1,n3)/3
!!Error: rdrg not defined
!!LLP user defined 
!         bustr(i)=rdrg*cff2
!         bvstr(i)=rdrg*cff1
!       END DO
!!#endif
       endif

       END SUBROUTINE set_vbc
