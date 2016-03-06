!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEEP_WATER(IP, ACLOC, IMATRA, IMATDA)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)      :: IP

         REAL(rkind), INTENT(OUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind), INTENT(IN)  :: ACLOC(MSC,MDC)

         IMATRA = 0.d0
         IMATDA = 0.d0

         IF (MESIN == 1) THEN
           WRITE(*,*) 'DOING ST4'
           CALL ST4_PRE(IP, ACLOC, IMATRA, IMATDA)
         ELSE IF (MESIN == 2) THEN
           WRITE(*,*) 'DOING WAM'
           CALL ECMWF_PRE(IP, ACLOC, IMATRA, IMATDA)
         ELSE IF (MESIN == 3) THEN
           WRITE(*,*) 'DOING SWAN SHIT'
           CALL CYCLE3_PRE(IP, ACLOC, IMATRA, IMATDA)
         ENDIF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************

