#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ECMWF_PRE (IP, ACLOC, IMATRA, IMATDA, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IP

         REAL(rkind), INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: SSINE(MSC,MDC),DSSINE(MSC,MDC), SSINL(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: SSDS(MSC,MDC),DSSDS(MSC,MDC)
         REAL(rkind), INTENT(OUT)   :: SSNL4(MDC,MSC),DSSNL4(MDC,MSC)

         INTEGER      :: IS, ID

         REAL(rkind)  :: VEC2RAD
         REAL(rkind)  :: WIND10
         REAL(rkind)  :: FPM,WINDTH,TEMP
         REAL(rkind)  :: SC, SP, JAC

         REAL(rkind)                :: FL3(MDC,MSC), FL(MDC,MSC), SL(MDC,MSC)

         SSINE = ZERO; DSSINE = ZERO
         SSDS  = ZERO; DSSDS  = ZERO
         SSNL4 = ZERO; DSSNL4 = ZERO
         SSINL = ZERO

         DO IS = 1, MSC
           DO ID = 1, MDC
             FL3(ID,IS) = ACLOC(IS,ID) * PI2 * SPSIG(IS)
           END DO
         END DO

         THWOLD(IP) = THWNEW(IP)
         U10NEW(IP) = MAX(TWO,SQRT(WINDXY(IP,1)**2+WINDXY(IP,2)**2))*WINDFAC
         Z0NEW(IP) = Z0OLD(IP)
         THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))

         CALL WAM_PRE (IP, FL3, FL, SL, SSDS, DSSDS, SSNL4, DSSNL4, SSINE, DSSINE)

         DO ID = 1, MDC
           DO IS = 1, MSC 
             JAC = ONE/PI2/SPSIG(IS)
             IMATRA(IS,ID) = SL(ID,IS)*JAC
             IMATDA(IS,ID) = FL(ID,IS)/PI2
           ENDDO
         ENDDO

         IF (.NOT. LINID) THEN
           CALL SET_WIND( IP, WIND10, WINDTH )
           CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
           CALL SIN_LIN_CAV(IP,WINDTH,FPM,SSINL)
           IMATRA = IMATRA + SSINL
         ENDIF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ECMWF_POST(IP,ACLOC)
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)           :: IP
         REAL(rkind), INTENT(INOUT)    :: ACLOC(MSC,MDC)
         INTEGER                       :: IS, ID
         REAL(rkind)                   :: VEC2RAD, FPM
         REAL(rkind)                   :: IMATRA(MSC,MDC)
         REAL(rkind)                :: FL3(MDC,MSC), FL(MDC,MSC), SL(MDC,MSC)
         THWOLD(IP) = THWNEW(IP)
         THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))
         U10NEW(IP) = MAX(TWO,SQRT(WINDXY(IP,1)**2+WINDXY(IP,2)**2)) * WINDFAC
         Z0NEW(IP) = Z0OLD(IP)
         DO IS = 1, MSC
           DO ID = 1, MDC
             FL3(ID,IS) =  ACLOC(IS,ID) * PI2 * SPSIG(IS)
           END DO
         END DO
         CALL WAM_POST (IP, FL3)
         DO IS = 1, MSC
           DO ID = 1, MDC
             ACLOC(IS,ID) = FL3(ID,IS) / PI2 / SPSIG(IS)
           END DO
         END DO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
