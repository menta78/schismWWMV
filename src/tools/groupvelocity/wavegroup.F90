!**********************************************************************
!*                                                                    *
!**********************************************************************
      PROGRAM wavegroup
         IMPLICIT NONE
         REAL(8) :: DEPLOC, SIGIN
         REAL(8) :: WVC, WVK, WVCG2, WN
         REAL(8) :: SGDLS , AUX1, AUX2
         REAL(8) :: WKDEP
! 
         WVC = 0.
         WVK = 0.
         WVCG2 = 0.
         WN = 0.

         WRITE(*,*) 'PLEASE PROVIDE SIG and DEP'
         READ(*,*) sigin, deploc

         IF (SIGIN .LT. 0.00000001) THEN
            WN = 0.
            WVK=10.
            WVCG2=0.
            RETURN
         END IF

         IF (DEPLOC .GT. 0.01) THEN
            SGDLS = SIGIN*SIGIN*DEPLOC/9.81
            AUX1 = 1.0+0.6522*SGDLS+0.4622*(SGDLS**2.0)+0.0864*(SGDLS**4.0)+0.0675*(SGDLS**5.0)
            AUX2 = 1.0/(SGDLS+1.0/AUX1)
            WVC = SQRT(AUX2*9.81*DEPLOC)
            WVK = SIGIN/WVC
            WKDEP = WVK*DEPLOC
            IF (WKDEP > 13.0) THEN
               WN = 0.5
            ELSE
               WN = 0.5*(1.+2.*WKDEP/SINH(MIN(2.*200.,2.*WKDEP)))
            END IF
            WVCG2 = WN*WVC
          ELSE
            WVC  = 0.
            WVK  = 10.
            WVCG2 = 0.
          END IF

          write(*,*) wvc, wvk, wvcg2
      END PROGRAM
!**********************************************************************
!*                                                                    *
!**********************************************************************

