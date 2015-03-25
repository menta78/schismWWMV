      SUBROUTINE sed_avalanching(dhnd)
!--------------------------------------------------------------------!
! This routine updates the bathymetry to account for avalanching     !
! The bed slopes are computed at each node. When a critical slope for!
! element is exceeded, bathymetry of each element node is modified   !
! in order to obtain a slope lower than critical threshold. The      !
! method used here conserves the volume.                             !
! Adapted from filter.f (SAND2D, A. Fortunato)                       !
!                                                                    !
! Currently it does not take account for modification bed sediment   !
! characteristics related to modification of bathymetry              !
!                                                                    !
! Author: florian ganthy (fganthy@lnec.pt ; florian.ganthy@gmail.com)!
! Date:   2013/01/16                                                 !
!                                                                    !
! History: 2013/01 - F.Ganthy : Modification of wet/dry element      !
!                    consideration to be more physical.              !
!                                                                    !
!--------------------------------------------------------------------!

      USE elfe_glbl, ONLY: rkind,dl,dp,nm,np,npa,nea,idry,errmsg,    &
                           area,xnd,ynd
      USE elfe_msgp, ONLY: myrank,parallel_abort,exchange_p2d,       &
                           exchange_e2d
      USE sed_mod,   ONLY: dry_slope_cr,wet_slope_cr

      IMPLICIT NONE
      SAVE

!- Local variables --------------------------------------------------!

      REAL(rkind), DIMENSION(npa), INTENT(inout) :: dhnd
      

      INTEGER :: i,j,iter,iflag,n1,n2,n3,ndry
      INTEGER,PARAMETER :: maxiter = 25 !Maximum number of iterations

      REAL(rkind) :: epsi,slope_cr,slope,h1,h2,h3,vec1,vec2,vec3,    &
                     m11,m12,m13,m21,m22,m23,m31,m32,m33,det,        &
                     h1p,h2p,h3p

      REAL(rkind), DIMENSION(:)   :: dp0(npa),dp1(npa),area2(npa),   &
                                     dph(npa),tmp(npa)
      REAL(rkind), DIMENSION(:,:) :: dpdxy_el(nea,2)
      

!- Start Statement --------------------------------------------------!
      IF(myrank.EQ.0) WRITE(16,*)'Entering sed_avalanching'


!--------------------------------------------------------------------!
! * Compute bed changes due to sediment transport
!--------------------------------------------------------------------!
      DO i=1,np
        dp0(i) = dp(i)+dhnd(i)
      ENDDO
      dp1 = dp0
      CALL exchange_p2d(dp0)
      CALL exchange_p2d(dp1)

!--------------------------------------------------------------------!
! * Compute node area
!--------------------------------------------------------------------!
      area2 = 0.0d0
      tmp = 0.0d0
      DO i=1,nea
        DO j=1,3
          area2(nm(i,j)) = area2(nm(i,j)) + (1.0d0/3.0d0)*area(i)
          tmp(nm(i,j)) = tmp(nm(i,j)) + 1.0d0
        ENDDO ! End loop j
      ENDDO ! End loop nea
      CALL exchange_p2d(area2)
!--------------------------------------------------------------------!
! * Start iterative procedure
!--------------------------------------------------------------------!
      iter  = 0
      iflag = 1
      epsi  = 0.01d0
      DO WHILE (iflag.EQ.1)
        iflag = 0
        iter  = iter+1
        dpdxy_el = 0.0d0
        DO i=1,nea
!--------------------------------------------------------------------!
! * Compute bed slope at element center
!--------------------------------------------------------------------!
          DO j=1,3
            dpdxy_el(i,1) = dpdxy_el(i,1)+dp1(nm(i,j))*dl(i,j,1)
            dpdxy_el(i,2) = dpdxy_el(i,2)+dp1(nm(i,j))*dl(i,j,2)
          ENDDO
          slope = dsqrt(dpdxy_el(i,1)*dpdxy_el(i,1)+                 &
          &             dpdxy_el(i,2)*dpdxy_el(i,2))
!--------------------------------------------------------------------!
! * Testing for critical slope value (dry or wet)
!   Here we consider that element is dry only if its three nodes are 
!   dry. In the former version, the test was done over idry_e (=1), 
!   but this potentially overestimated the slumping at the wet-dry
!   limit
!--------------------------------------------------------------------!
          ndry = idry(nm(i,1))+idry(nm(i,3))+idry(nm(i,3))
          IF(ndry.EQ.3) THEN
            slope_cr = dry_slope_cr
          ELSE
            slope_cr = wet_slope_cr
          ENDIF

          IF(slope-slope_cr.LE.epsi) CYCLE
!--------------------------------------------------------------------!
! * Preparation of system equation
!--------------------------------------------------------------------!
          n1 = nm(i,1)
          n2 = nm(i,2)
          n3 = nm(i,3)
          h1 = dp1(n1)
          h2 = dp1(n2)
          h3 = dp1(n3)
          m11 = area2(n1)
          m12 = area2(n2)
          m13 = area2(n3)
          m21 = ynd(n2)-ynd(n3)
          m22 = ynd(n3)-ynd(n1)
          m23 = ynd(n1)-ynd(n2)
          m31 = xnd(n3)-xnd(n2)
          m32 = xnd(n1)-xnd(n3)
          m33 = xnd(n2)-xnd(n1)
          vec1 = h1*m11 + h2*m12 + h3*m13
          vec2 = slope_cr/slope * (h1*m21 + h2*m22 + h3*m23)
          vec3 = slope_cr/slope * (h1*m31 + h2*m32 + h3*m33)
!--------------------------------------------------------------------!
! * Solving the system by Cramer's rule
!--------------------------------------------------------------------!
          det = m11*(m22*m33-m32*m23)-                               &
          &     m12*(m21*m33-m31*m23)+                               &
          &     m13*(m21*m32-m31*m22)

          IF(det.EQ.0.0d0) THEN
            WRITE(errmsg,*)'SED_AVALANCHING: det=0.0'
            CALL parallel_abort(errmsg)
          ENDIF
!--------------------------------------------------------------------!
! * Compute new depth at nodes
!--------------------------------------------------------------------!
          h1p = (vec1*(m22*m33-m32*m23)-                             &
          &     m12*(vec2*m33-vec3*m23)+                             &
          &     m13*(vec2*m32-vec3*m22))/det

          h2p = (m11*(vec2*m33-vec3*m23)-                            &
          &     vec1*(m21*m33-m31*m23)+                              &
          &     m13*(m21*vec3-m31*vec2))/det

          h3p = (m11*(m22*vec3-m32*vec2)-                            &
          &     m12*(m21*vec3-m31*vec2)+                             &
          &     vec1*(m21*m32-m31*m22))/det
!--------------------------------------------------------------------!
! * Apply new depth to temporary bathymetry
!--------------------------------------------------------------------!
          dp1(n1) = h1p
          dp1(n2) = h2p
          dp1(n3) = h3p
          
          iflag = 1
          IF (iter.GE.maxiter) iflag=0

        ENDDO ! End loop on nea
      ENDDO ! End loop on iflag

!--------------------------------------------------------------------!
! * Apply depth changes
!--------------------------------------------------------------------!
      DO i=1,np
        dhnd(i) = dhnd(i)+(dp1(i)-dp0(i))
        IF(isnan(dhnd(i)).EQV..TRUE.) THEN
          WRITE(errmsg,*) 'Avalanching: dhnd is NaN',myrank,dhnd(i), &
                          dp1(i),dp0(i)
          CALL parallel_abort(errmsg)
        ENDIF
      ENDDO
     CALL exchange_p2d(dhnd)
!--------------------------------------------------------------------!
! * Write number of iteration within mirror.out
!--------------------------------------------------------------------!
      IF(myrank.EQ.0) THEN
        WRITE(16,*)'Number of avalanching iterations:',iter
      ENDIF

!--------------------------------------------------------------------!
      IF(myrank.EQ.0) WRITE(16,*)'Leaving sed_avalanching'
!--------------------------------------------------------------------!
      END SUBROUTINE sed_avalanching
