!**********************************************************************
!*   This is source code of MAKIN formulation                         *
!*   It is not incorporated in the source code of WWM but could       *
!*   be put back. The ICOMP stuff would then needs to be cleaned up   *
!*   and probably other stuff.                                        *
!**********************************************************************
      SUBROUTINE SIN_MAKIN(IP, WIND10, WINDTH, KMESPC, ETOT, WALOC, PHI, DPHIDN, SSINE)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP
         REAL(rkind)   , INTENT(IN)    :: WIND10, WINDTH
         REAL(rkind)   , INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind)   , INTENT(OUT)   :: SSINE(NUMSIG,NUMDIR)
         REAL(rkind)   , INTENT(INOUT) :: PHI(NUMSIG,NUMDIR), DPHIDN(NUMSIG,NUMDIR)

         INTEGER             :: IS, ID
         REAL(rkind)                :: AUX1, AUX2
         REAL(rkind)                :: SWINB, CINV, SIGMA
         REAL(rkind)                :: COSWW, THETA
         REAL(rkind)                :: NC_MK, MC_MK, MBETA, RMK
         REAL(rkind)                :: KMESPC, ETOT, DS, ELOCAL
         REAL(rkind)                :: STEEPLOCAL
         REAL(rkind)                :: ALOCAL, CPHASE, CYS, ATTC

         LOGICAL             :: LATT, LOPP
!
! PARAMETER FROM MAKIN 1999
!
         MC_MK = 0.3_rkind
         NC_MK = 5.0_rkind
         LOPP  = .FALSE.
         CYS   = -25._rkind  ! Opposing wind attenuation.
         LATT  = .FALSE.
         ATTC  = -10.0_rkind  ! Attenuation coefficient
         MBETA =  32.0_rkind   ! See orignial Paper A GU OCEANS VOL. 104, No.: C4, April 1999 and see Makin & Stam 2003 (KNMI)

         DO IS = 1, NUMSIG
           CINV =  WK(IS,IP) / SPSIG(IS) 
           SIGMA = SPSIG(IS)
           IF (WIND10 .LE. THR) THEN
             AUX1 = 0.0_rkind
           ELSE
             AUX1  = 1._rkind/CINV/WIND10
           END IF
           AUX2  = UFRIC(IP)*CINV
           RMK = 1 - MC_MK * AUX1 ** NC_MK
           DO ID = 1, NUMDIR
             THETA  = SPDIR(ID)
             COSWW  = MyCOS(THETA-WINDTH)
             IF (LATT) THEN
               IF (RMK .GE. 0.0_rkind) THEN
                 SWINB = MBETA * RMK * RHOAW * AUX2**2 * COSWW * ABS(COSWW) * SIGMA
               END IF
               IF (COSWW * ABS(COSWW) .GE. 0.0_rkind .AND. RMK .LT. 0.0_rkind) THEN
                 SWINB = MAX(ATTC,MBETA*RMK) * RHOAW *  AUX2**2 * COSWW * ABS(COSWW) * SIGMA
               END IF
             ELSE
               IF (RMK .GT. ZERO) THEN
                 SWINB = MBETA * RMK * RHOAW * AUX2**2 * COSWW * ABS(COSWW) * SIGMA
               ELSE
                 SWINB = ZERO
               END IF
             END IF
             IF (COSWW * ABS(COSWW) .LE. ZERO) THEN
               IF (LOPP) THEN
                 CPHASE      = 0.0_rkind/CINV
                 IF (IS .EQ. 1) DS = SPSIG(IS)
                 IF (IS .GT. 1) DS = SPSIG(IS) - SPSIG(IS-1)
                 IF (IS .EQ. 1) ELOCAL = ONEHALF * WALOC(IS,ID) * SPSIG(IS) * SPSIG(IS) ! Simpson
                 IF (IS .GT. 1) ELOCAL = ONEHALF * ( WALOC(IS,ID) * SPSIG(IS) + WALOC(IS-1,ID) * SPSIG(IS-1) ) * DS 
                 ALOCAL      = SQRT(8.0_rkind*ELOCAL)
                 STEEPLOCAL  = ALOCAL  * WK(IS,IP)
                 SWINB       = CYS * RHOAW * STEEPLOCAL * STEEPLOCAL * (ONE - ((WIND10 * COSWW)/CPHASE) ) **2 * SIGMA
               ELSE
                 SWINB = ZERO
               END IF
             END IF

             SSINE(IS,ID)   = SWINB * WALOC(IS,ID)

             IF (ICOMP .GE. 2) THEN
               IF (SWINB .LT. 0) THEN
                 DPHIDN(IS,ID) = - SWINB
               ELSE
                 PHI(IS,ID) =  PHI(IS,ID) + SSINE(IS,ID)
               END IF
             ELSE IF (ICOMP .LT. 2) THEN
               IF (SWINB .LT. 0) THEN
                 DPHIDN(IS,ID) = SWINB
                 PHI(IS,ID) = PHI(IS,ID) + SSINE(IS,ID)
               ELSE
                 PHI(IS,ID) = PHI(IS,ID) + SSINE(IS,ID)
               END IF
             END IF

           END DO
         END DO

         RETURN
      END SUBROUTINE
