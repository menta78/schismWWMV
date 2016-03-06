!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEEP_WATER(IP, ACLOC, IMATRA, IMATDA, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)      :: IP
         REAL(rkind), INTENT(IN)  :: ACLOC(MSC,MDC)

         REAL(rkind), INTENT(OUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind), INTENT(OUT) :: SSINL(MSC,MDC)
         REAL(rkind), INTENT(OUT) :: SSINE(MSC,MDC),DSSINE(MSC,MDC)
         REAL(rkind), INTENT(OUT) :: SSNL4(MSC,MDC),DSSNL4(MSC,MDC)
         REAL(rkind), INTENT(OUT) :: SSDS(MSC,MDC),DSSDS(MSC,MDC)

         IF (ISOURCE == 1) THEN
           WRITE(*,*) 'DOING ST4'
           CALL ST4_PRE(IP, ACLOC, IMATRA, IMATDA, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
         ELSE IF (ISOURCE == 2) THEN
           WRITE(*,*) 'DOING WAM'
           CALL ECMWF_PRE(IP, ACLOC, IMATRA, IMATDA, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
           WRITE(*,*) SUM(IMATRA), SUM(IMATDA), 'DEEP WATER' 
         ELSE IF (ISOURCE == 3) THEN
           WRITE(*,*) 'DOING SWAN SHIT'
           CALL CYCLE3_PRE(IP, ACLOC, IMATRA, IMATDA, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
         ENDIF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
