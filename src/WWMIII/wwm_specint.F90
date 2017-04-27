#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SEMI_IMPLICIT_INTEGRATION(IP,DT,WALOC)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IP
      REAL(rkind), INTENT(IN) :: DT
      REAL(rkind), INTENT(INOUT) :: WALOC(NUMSIG,NUMDIR)
      INTEGER       :: IS, ID
      REAL(rkind)   :: NEWDAC, SSBR(NUMSIG,NUMDIR)
      REAL(rkind)   :: PHI(NUMSIG,NUMDIR), DPHIDN(NUMSIG,NUMDIR), ACOLD(NUMSIG,NUMDIR)

      ACOLD = WALOC
      CALL COMPUTE_PHI_DPHI(IP,WALOC,PHI,DPHIDN)
      DO IS = 1, NUMSIG
        DO ID = 1, NUMDIR
          NEWDAC = PHI(IS,ID) * DT / (ONE-DT*MIN(ZERO,DPHIDN(IS,ID)))
          WALOC(IS,ID) = MAX( ZERO, ACOLD(IS,ID) + NEWDAC )
        END DO
      END DO
      !CALL POST_INTEGRATION(IP,WALOC)
      IF (LLIMT) CALL LIMITER(IP,ACOLD,WALOC)
      IF (LMAXETOT) CALL BREAK_LIMIT(IP,WALOC,SSBR)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_PHI_DPHI(IP,WALOC,PHI,DPHIDN)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN)      :: IP
      REAL(rkind), INTENT(IN)  :: WALOC(NUMSIG,NUMDIR)
      REAL(rkind), INTENT(OUT) :: PHI(NUMSIG,NUMDIR), DPHIDN(NUMSIG,NUMDIR)
      REAL(rkind)   :: SSINL(NUMSIG,NUMDIR)
      REAL(rkind)   :: SSINE(NUMSIG,NUMDIR),DSSINE(NUMSIG,NUMDIR)
      REAL(rkind)   :: SSDS(NUMSIG,NUMDIR),DSSDS(NUMSIG,NUMDIR)
      REAL(rkind)   :: SSNL3(NUMSIG,NUMDIR),DSSNL3(NUMSIG,NUMDIR)
      REAL(rkind)   :: SSNL4(NUMSIG,NUMDIR),DSSNL4(NUMSIG,NUMDIR)
      REAL(rkind)   :: SSBR(NUMSIG,NUMDIR),DSSBR(NUMSIG,NUMDIR)
      REAL(rkind)   :: SSBRL(NUMSIG,NUMDIR)
      REAL(rkind)   :: SSBF(NUMSIG,NUMDIR),DSSBF(NUMSIG,NUMDIR)
      REAL(rkind)   :: HS,TM01,TM02,TM10,KLM,WLM
#ifdef DEBUG
      REAL(rkind)   :: SSINL_WW3(NUMSIG,NUMDIR), SSINE_WW3(NUMSIG,NUMDIR)
      REAL(rkind)   :: SSBRL_WW3(NUMSIG,NUMDIR), SSDS_WW3(NUMSIG,NUMDIR)
      REAL(rkind)   :: SSNL4_WW3(NUMSIG,NUMDIR), SSBF_WW3(NUMSIG,NUMDIR)
      REAL(rkind)   :: SSBR_WW3(NUMSIG,NUMDIR)
#endif

      DPHIDN = ZERO; PHI    = ZERO
      SSINE  = ZERO; DSSINE = ZERO
      SSDS   = ZERO; DSSDS  = ZERO
      SSNL3  = ZERO; DSSNL3 = ZERO
      SSNL4  = ZERO; DSSNL4 = ZERO
      SSBR   = ZERO; DSSBR  = ZERO
      SSBF   = ZERO; DSSBF  = ZERO
      SSBRL  = ZERO
      SSINL  = ZERO

      IF (IOBP(IP) .NE. 0 .AND. .NOT. LSOUBOUND) THEN
        RETURN
      ELSE IF (LSOUBOUND .AND. IOBP(IP) .EQ. 2) THEN
        RETURN
      ENDIF

#ifdef DEBUG
      IF (IP .eq. TESTNODE) THEN
         WRITE(740+myrank,*) 'Before integration'
         WRITE(740+myrank,*) 'sum(WALOC)=', sum(WALOC)
         CALL MEAN_PARAMETER(IP,WALOC,NUMSIG,HS,TM01,TM02,TM10,KLM,WLM)
         WRITE(740+myrank,*) 'HS=', HS, ' TM01=', TM01
      END IF
#endif
!
      CALL DEEP_WATER(IP, WALOC, PHI, DPHIDN, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
!
#ifdef DEBUG
      IF (IP .eq. TESTNODE) THEN
         WRITE(740+myrank,*) 'sum(SSINE)=', sum(SSINE), ' sum(DSSINE)=', sum(DSSINE)
      END IF
#endif
!
      IF (ISHALLOW(IP) .EQ. 1) CALL SHALLOW_WATER(IP, WALOC, PHI, DPHIDN, SSBR, DSSBR, SSBF, DSSBF, SSBRL, SSNL3, DSSNL3)
!
#ifdef DEBUG
      IF (IP .eq. TESTNODE) THEN
         WRITE(740+myrank,*) 'After integration ISOURCE=', ISOURCE
         WRITE(740+myrank,*) 'sum(WALOC)=', sum(WALOC)
         CALL MEAN_PARAMETER(IP,WALOC,NUMSIG,HS,TM01,TM02,TM10,KLM,WLM)
         WRITE(740+myrank,*) 'HS=', HS, ' TM01=', TM01
         CALL CONVERT_VS_WWM_TO_WW3(IP, SSINL, SSINL_WW3)
         WRITE(740+myrank,*) 'WW3 : LINEAR INPUT =', SUM(SSINL_WW3), MINVAL(SSINL_WW3), MAXVAL(SSINL_WW3)
         CALL CONVERT_VS_WWM_TO_WW3(IP, SSINE, SSINE_WW3)
         WRITE(740+myrank,*) 'WW3 : EXP. INPUT   =', SUM(SSINE_WW3), MINVAL(SSINE_WW3), MAXVAL(SSINE_WW3)
         CALL CONVERT_VS_WWM_TO_WW3(IP, SSBRL, SSBRL_WW3)
         WRITE(740+myrank,*) 'WW3 : BREAK LIMIT  =', SUM(SSBRL_WW3), MINVAL(SSBRL_WW3), MAXVAL(SSBRL_WW3)
         CALL CONVERT_VS_WWM_TO_WW3(IP, SSDS, SSDS_WW3)
         WRITE(740+myrank,*) 'WW3 : WHITECAP     =', SUM(SSDS_WW3), MINVAL(SSDS_WW3), MAXVAL(SSDS_WW3)
         CALL CONVERT_VS_WWM_TO_WW3(IP, SSNL4, SSNL4_WW3)
         WRITE(740+myrank,*) 'WW3 : SNL4         =', SUM(SSNL4_WW3), MINVAL(SSNL4_WW3), MAXVAL(SSNL4_WW3)
         CALL CONVERT_VS_WWM_TO_WW3(IP, SSBF, SSBF_WW3)
         WRITE(740+myrank,*) 'WW3 : BOT. FRIC.   =', SUM(SSBF_WW3), MINVAL(SSBF_WW3), MAXVAL(SSBF_WW3)
         CALL CONVERT_VS_WWM_TO_WW3(IP, SSBR, SSBR_WW3)
         WRITE(740+myrank,*) 'WW3 : BREAKING     =', SUM(SSBR_WW3), MINVAL(SSBR_WW3), MAXVAL(SSBR_WW3)
         
         WRITE(740+myrank,*) 'WAVE ACTION      =', SUM(WALOC), MINVAL(WALOC), MAXVAL(WALOC)
         WRITE(740+myrank,*) 'LINEAR INPUT     =', SUM(SSINL), MINVAL(SSINL), MAXVAL(SSINL)
         WRITE(740+myrank,*) 'BREAKING LIMITER =', SUM(SSBRL), MINVAL(SSBRL), MAXVAL(SSBRL)
         WRITE(740+myrank,*) 'EXP INPUT(SSINE) =', SUM(SSINE),  MINVAL(SSINE),  MAXVAL(SSINE)
         WRITE(740+myrank,*) 'EXP INPUT(DSSINE)=', SUM(DSSINE), MINVAL(DSSINE), MAXVAL(DSSINE)
         WRITE(740+myrank,*) 'WHITECAP(SSDS)   =', SUM(SSDS),  MINVAL(SSDS),  MAXVAL(SSDS)
         WRITE(740+myrank,*) 'WHITECAP(DSSDS)  =', SUM(DSSDS), MINVAL(DSSDS), MAXVAL(DSSDS)
         WRITE(740+myrank,*) 'SNL4(SSNL4)      =', SUM(SSNL4), MINVAL(SSNL4), MAXVAL(SSNL4)
         WRITE(740+myrank,*) 'SNL4(DSSNL4)     =', SUM(DSSNL4), MINVAL(DSSNL4), MAXVAL(DSSNL4)
         WRITE(740+myrank,*) 'BOT. FRIC(SSBF)  =', SUM(SSBF),  MINVAL(SSBF),  MAXVAL(SSBF)
         WRITE(740+myrank,*) 'BOT. FRIC(DSSBF) =', SUM(DSSBF), MINVAL(DSSBF), MAXVAL(DSSBF)
         WRITE(740+myrank,*) 'BREAKING(SSBR)   =', SUM(SSBR),  MINVAL(SSBR),  MAXVAL(SSBR)
         WRITE(740+myrank,*) 'BREAKING(DSSBR)  =', SUM(DSSBR), MINVAL(DSSBR), MAXVAL(DSSBR)
         WRITE(740+myrank,*) 'TOTAL SRC(PHI)=', SUM(PHI), MINVAL(PHI), MAXVAL(PHI)
         WRITE(740+myrank,*) 'TOTAL SRC(DPHIDN)=', SUM(DPHIDN), MINVAL(DPHIDN), MAXVAL(DPHIDN)
      END IF
#endif

#ifdef DEBUG_SOURCE_TERM
      WRITE(*,'(A20,6E20.10)') 'WAVE ACTION', SUM(WALOC), MINVAL(WALOC), MAXVAL(WALOC)
      WRITE(*,'(A20,6E20.10)') 'LINEAR INPUT', SUM(SSINL), MINVAL(SSINL), MAXVAL(SSINL)
      WRITE(*,'(A20,6E20.10)') 'BREAKING LIMITER', SUM(SSBRL), MINVAL(SSBRL), MAXVAL(SSBRL)
      WRITE(*,'(A20,6E20.10)') 'EXP INPUT', SUM(SSINE), SUM(DSSINE), MINVAL(SSINE), MAXVAL(SSINE), MINVAL(DSSINE), MAXVAL(DSSINE)
      WRITE(*,'(A20,6E20.10)') 'WHITECAP', SUM(SSDS), SUM(DSSDS), MINVAL(SSDS), MAXVAL(SSDS), MINVAL(DSSDS), MAXVAL(DSSDS)
      WRITE(*,'(A20,6E20.10)') 'SNL4', SUM(SSNL4), SUM(DSSNL4), MINVAL(SSNL4), MAXVAL(SSNL4), MINVAL(DSSNL4), MAXVAL(DSSNL4)
      WRITE(*,'(A20,6E20.10)') 'BOTTOM FRICTION', SUM(SSBF), SUM(DSSBF), MINVAL(SSBF), MAXVAL(SSBF), MINVAL(DSSBF), MAXVAL(DSSBF)
      WRITE(*,'(A20,6E20.10)') 'BREAKING', SUM(SSBR), SUM(DSSBR), MINVAL(SSBR), MAXVAL(SSBR), MINVAL(DSSBR), MAXVAL(DSSBR)
      WRITE(*,'(A20,6E20.10)') 'TOTAL SOURCE TERMS', SUM(PHI), SUM(DPHIDN), MINVAL(PHI), MAXVAL(PHI), MINVAL(DPHIDN), MAXVAL(DPHIDN)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CONVERT_VS_WWM_TO_WW3(IP, VS_WWM, VS_WW3)
      USE DATAPOOL
      IMPLICIT NONE
      integer IP
      REAL(rkind), intent(in) :: VS_WWM(NUMSIG,NUMDIR)
      REAL(rkind), intent(out) :: VS_WW3(NUMSIG,NUMDIR)
      INTEGER ID,IS
      DO ID=1,NUMDIR
        DO IS=1,NUMSIG
          VS_WW3(IS,ID) = CG(IS,IP) * VS_WWM(IS,ID)
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE POST_INTEGRATION(IP, WALOC)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)       :: IP
         REAL(rkind),INTENT(INOUT) :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind)               :: SSINE(NUMSIG,NUMDIR),DSSINE(NUMSIG,NUMDIR), SSINL(NUMSIG,NUMDIR)
         REAL(rkind)               :: SSDS(NUMSIG,NUMDIR),DSSDS(NUMSIG,NUMDIR)

         IF (ISOURCE == 1) THEN
           CALL ST4_POST(IP, WALOC, SSINE, DSSINE, SSDS, DSSDS, SSINL)
         ELSE IF (ISOURCE == 2) THEN
           CALL ECMWF_POST(IP, WALOC)
         ELSE IF (ISOURCE == 3) THEN
!2do write some post code for cycle3
         ENDIF
      END SUBROUTINE POST_INTEGRATION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SOURCES_EXPLICIT
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP
         REAL(rkind) :: WALOC(NUMSIG,NUMDIR)

         IF (SMETHOD .gt. 0) THEN
           DO IP = 1, MNP
             IF (IOBDP(IP) .EQ. 1) THEN
               WALOC = AC2(:,:,IP)
               CALL SEMI_IMPLICIT_INTEGRATION(IP,DT4S,WALOC)
               AC2(:,:,IP) = WALOC
             ENDIF
           ENDDO
         ENDIF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SOURCES_IMPLICIT
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER     :: IP
         REAL(rkind) :: WALOC(NUMSIG,NUMDIR),PHI(NUMSIG,NUMDIR),DPHIDN(NUMSIG,NUMDIR)

         DO IP = 1, MNP
           WALOC = AC2(:,:,IP)
           CALL COMPUTE_PHI_DPHI(IP,WALOC,PHI,DPHIDN)
           PHIA(:,:,IP)    = PHI
           DPHIDNA(:,:,IP) = DPHIDN
         ENDDO

#ifdef DEBUG_SOURCE_TERM
         WRITE(*,*) 'SOURCES_IMPLICIT', SUM(PHIA), SUM(DPHIDNA)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
         SUBROUTINE LIMITER(IP,ACOLD,WALOC)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IP
         INTEGER                    :: IS, ID
         REAL(rkind), INTENT(INOUT) :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(IN)    :: ACOLD(NUMSIG,NUMDIR)
         REAL(rkind)                :: NEWDAC, OLDAC, NEWAC, DELT, XIMP, DELFL(NUMSIG)
         REAL(rkind)                :: MAXDAC, CONST, SND, DELT5, USFM, eFric, PHILMAXDAC

!         CONST  = PI2**2*3.0*1.0E-7*DT4S*SPSIG(NUMSIG)
!         SND    = PI2*5.6*1.0E-3
         DELT   = DT4S
!         XIMP   = 1._rkind
!         DELT5  = XIMP*DELT
         DELFL  = COFRM4*DELT
 
         DO IS = 1, NUMSIG
           PHILMAXDAC = 0.0081*LIMFAK/(TWO*SPSIG(IS)*WK(IS,IP)**3*CG(IS,IP))
           IF (ISOURCE .EQ. 1) THEN
              USFM   = UFRIC(IP)*MAX(FMEANWS(IP),FMEAN(IP))
              MAXDAC = MAX(PHILMAXDAC,USFM*DELFL(IS)/PI2/SPSIG(IS))
           ELSE IF (ISOURCE .EQ. 2) THEN
              USFM   = USNEW(IP)*MAX(FMEANWS(IP),FMEAN(IP))
              MAXDAC = MAX(PHILMAXDAC,USFM*DELFL(IS)/PI2/SPSIG(IS))
           ELSE IF (ISOURCE .EQ. 3) THEN
              MAXDAC = PHILMAXDAC
           END IF
           DO ID = 1, NUMDIR
             NEWAC  = WALOC(IS,ID)
             OLDAC  = ACOLD(IS,ID)
             NEWDAC = NEWAC - OLDAC
             NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)), NEWDAC)
             WALOC(IS,ID) = MAX( ZERO, OLDAC + NEWDAC )
           END DO
         END DO

         END SUBROUTINE LIMITER
!**********************************************************************
!*                                                                    *
!**********************************************************************
         SUBROUTINE BREAK_LIMIT(IP,WALOC,SSBRL)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)  :: IP

         REAL(rkind), INTENT(INOUT)  :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT)    :: SSBRL(NUMSIG,NUMDIR)

         REAL(rkind)                 :: EFTAIL, HS
         REAL(rkind)                 :: EMAX, RATIO, ETOT
         REAL(rkind)                 :: DINTSPEC

         ETOT   = 0.0
         EFTAIL = 1.0 / (TAIL_ARR(1)-1.0)

         ETOT = DINTSPEC(IP,WALOC)

         HS = 4.*SQRT(ETOT)

         !IF (LMONO_IN) HMAX(IP) = HMAX(IP) * SQRT(2.)

         EMAX = 1./16. * (HMAX(IP))**2

         IF (ETOT .GT. EMAX .AND. ETOT .GT. THR) THEN
           RATIO = EMAX/ETOT
           SSBRL = WALOC - RATIO * WALOC
           WALOC = RATIO * WALOC
         END IF
         END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
         SUBROUTINE RESCALE_SPECTRUM
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP, IS, ID
         REAL(rkind)    :: ETOTAL, EPOSIT
         REAL(rkind)    :: FACTOR
         LOGICAL :: ENEG

         DO IP = 1, MNP
            DO IS = 1, NUMSIG
               ETOTAL = ZERO 
               EPOSIT = ZERO 
               ENEG   = .FALSE.
               DO ID = 1, NUMDIR
                  ETOTAL = ETOTAL + AC2(IS,ID,IP)
                  IF (AC2(IS,ID,IP) > ZERO) THEN
                     EPOSIT = EPOSIT + AC2(IS,ID,IP)
                  ELSE
                     ENEG = .TRUE.
                  END IF
               END DO
               IF (ENEG) THEN
                  IF (EPOSIT .GT. VERYSMALL) THEN
                    FACTOR = ETOTAL/EPOSIT
                  ELSE 
                    FACTOR = 0.
                  END IF
                  DO ID = 1, NUMDIR
                     IF (AC2(IS,ID,IP) < ZERO) AC2(IS,ID,IP) = ZERO 
                     IF (FACTOR >= ZERO)  AC2(IS,ID,IP) = AC2(IS,ID,IP)*FACTOR
                     AC2(IS,ID,IP) = MAX(zero,AC2(IS,ID,IP))
                  END DO
               END IF
            END DO
         END DO
         END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
         SUBROUTINE SETSHALLOW
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER :: IP

         DO IP = 1, MNP
           IF (WK(1,IP)*DEP(IP) .LT. PI) THEN
             ISHALLOW(IP) = 1
           ELSE
             ISHALLOW(IP) = 0
           END IF
         END DO
         END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
         SUBROUTINE BREAK_LIMIT_ALL
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER              :: IP
         REAL(rkind)          :: HS
         REAL(rkind)          :: EMAX, RATIO, ETOT
         REAL(rkind)          :: DINTSPEC
         REAL(rkind)          :: ACLOC(NUMSIG, NUMDIR)
!      Print *, 'Passing BREAK_LIMIT_ALL'
         DO IP = 1, MNP
           ACLOC = AC2(:,:,IP)
           IF (ISHALLOW(IP) .EQ. 0) CYCLE
           ETOT = DINTSPEC(IP,ACLOC)
           HS = 4.*SQRT(ETOT)
           EMAX = 1./16. * (HMAX(IP))**2
!        WRITE(300,*) 'IP=', IP, ' HMAX=', HMAX(IP), ' DEP=', DEP(IP)
!        WRITE(300,*) '   ', IP, ' EMAX=', EMAX, ' ETOT=', ETOT
!        WRITE(300,*) '   ', IP, ' HS=', HS, ' BRHD=', BRHD

           IF (ETOT .GT. EMAX) THEN
             if(myrank==0) WRITE(300,*) '   break XP=', XP(IP)
             RATIO = EMAX/ETOT
             AC2(:,:,IP) = RATIO * ACLOC(:,:)
             AC1(:,:,IP) = RATIO * ACLOC(:,:)
           END IF
         END DO
         END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************

