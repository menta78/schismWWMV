!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEEP_WATER(IP, WALOC, PHI, DPHIDN, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)      :: IP
         REAL(rkind), INTENT(IN)  :: WALOC(NUMSIG,NUMDIR)

         REAL(rkind), INTENT(OUT) :: PHI(NUMSIG,NUMDIR), DPHIDN(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT) :: SSINL(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT) :: SSINE(NUMSIG,NUMDIR),DSSINE(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT) :: SSNL4(NUMSIG,NUMDIR),DSSNL4(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT) :: SSDS(NUMSIG,NUMDIR),DSSDS(NUMSIG,NUMDIR)

         IF (ISOURCE == 1) THEN
#ifdef DEBUG_SOURCE_TERM
           WRITE(STAT%FHNDL,*) 'DOING ST4'
#endif
           CALL ST4_PRE(IP, WALOC, PHI, DPHIDN, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
         ELSE IF (ISOURCE == 2) THEN
#ifdef DEBUG_SOURCE_TERM
           WRITE(STAT%FHNDL,*) 'DOING WAM'
#endif
           CALL ECMWF_PRE(IP, WALOC, PHI, DPHIDN, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
#ifdef DEBUG_SOURCE_TERM
           WRITE(STAT%FHNDL,*) SUM(PHI), SUM(DPHIDN), 'DEEP WATER' 
#endif
         ELSE IF (ISOURCE == 3) THEN
#ifdef DEBUG_SOURCE_TERM
           WRITE(STAT%FHNDL,*) 'DOING CYCLE 3'
#endif
           CALL CYCLE3_PRE(IP, WALOC, PHI, DPHIDN, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
         ENDIF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
