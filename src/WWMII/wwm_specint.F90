#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
#define ACTIVATE_SMETHOD_5
#undef ACTIVATE_SMETHOD_5
      SUBROUTINE SOURCE_INT_EXP()

         USE DATAPOOL
#ifdef ST41
         USE W3SRC4MD_OLD, ONLY : LFIRSTSOURCE
#endif
#ifdef ST42
         USE W3SRC4MD, ONLY : LFIRSTSOURCE
#endif
         IMPLICIT NONE

         INTEGER        :: IP, IS, ID, I, K, M
         INTEGER        :: NIT_SIN, NIT_SDS, NIT_SNL4, NIT_SNL3, NIT_SBR, NIT_SBF, NIT_ALL
         INTEGER, SAVE  :: IFIRST 
         REAL(rkind)    :: ACLOC(MSC,MDC), IMATRA(MSC,MDC), IMATDA(MSC,MDC), SSBRL2(MSC,MDC)
         REAL(rkind)    :: DT4S_T, DT4S_E, DT4S_Q, DT4S_H, DT4S_TQ, DT4S_TS

         ISELECT = 10

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(IP,IS,ID,ACLOC)
         DO IP = 1, MNP
!           IF (IP_IS_STEADY(IP) .EQ. 1) CYCLE
           IF ((ABS(IOBP(IP)) .NE. 1 .AND. IOBP(IP) .NE. 3)) THEN
             IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
               ACLOC  = AC2(IP,:,:)
               IF (SMETHOD == 1) THEN
                 CALL INT_IP_STAT(IP, DT4S, LLIMT, ACLOC)
               ELSE IF (SMETHOD == 2) THEN
                 CALL RKS_SP3(IP, DT4S, LLIMT,ACLOC)
               ELSE IF (SMETHOD == 3) THEN
                 CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_DYN, NDYNITER, ACLOC, NIT_ALL)
               END IF
               !CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .TRUE.) ! Update everything based on the new spectrum ... recalc
               IF (LMAXETOT .AND. .NOT. LADVTEST .AND. ISHALLOW(IP) .EQ. 1) CALL BREAK_LIMIT(IP,ACLOC,SSBRL2)
               AC2(IP,:,:) = ACLOC
             ENDIF
           ELSE !Boundary node ... 
             IF (LSOUBOUND) THEN ! Source terms on boundary ...
               IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
                 ACLOC  = AC2(IP,:,:)
                 IF (SMETHOD == 1) THEN
                   CALL INT_IP_STAT(IP, DT4S, LLIMT, ACLOC)
                 ELSE IF (SMETHOD == 2) THEN
                   CALL RKS_SP3(IP, DT4S,LLIMT,ACLOC)
                 ELSE IF (SMETHOD == 3) THEN
                   CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_DYN, NDYNITER, ACLOC, NIT_ALL)
                 ENDIF
                 IF (SMETHOD .GT. 0) THEN
                   !CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .TRUE.) ! Update everything based on the new spectrum ... recalc
                   IF (LMAXETOT .AND. .NOT. LADVTEST .AND. ISHALLOW(IP) .EQ. 1) THEN
                     CALL BREAK_LIMIT(IP,ACLOC,SSBRL2)
                   ENDIF
                   AC2(IP,:,:) = ACLOC
                 ENDIF
               ENDIF
             ELSE
               ACLOC = AC2(IP,:,:)
               IF (LMAXETOT .AND. .NOT. LADVTEST .AND. ISHALLOW(IP) .EQ. 1) THEN ! limit wave height on the boundary ...
                 CALL BREAK_LIMIT(IP,ACLOC,SSBRL2)
                 AC2(IP,:,:) = ACLOC
               ENDIF
             ENDIF 
           ENDIF
           IF (LNANINFCHK) THEN 
             IF (SUM(ACLOC) .NE. SUM(ACLOC) ) THEN 
               WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, '   IN SOURCE TERM INTEGRATION'
               CALL WWM_ABORT('wwm_specint.F90 l.88')
             END IF
           ENDIF
           AC1(IP,:,:) = AC2(IP,:,:)
         ENDDO
#if defined ST41 || defined ST42
         LFIRSTSOURCE = .FALSE.
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SOURCE_INT_IMP_WWM

         USE DATAPOOL
#ifdef ST41
         USE W3SRC4MD_OLD, ONLY : LFIRSTSOURCE
#endif
#ifdef ST42
         USE W3SRC4MD, ONLY : LFIRSTSOURCE
#endif
         IMPLICIT NONE

         INTEGER :: IP

         REAL(rkind)    :: ACLOC(MSC,MDC)
         REAL(rkind)    :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

!$OMP WORKSHARE
         IMATDAA = 0.
         IMATRAA = 0.
!$OMP END WORKSHARE

         ISELECT = 10

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(IP,ACLOC,IMATDA,IMATRA)
         DO IP = 1, MNP
           IF ((ABS(IOBP(IP)) .NE. 1 .AND. IOBP(IP) .NE. 3)) THEN
             IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
               ACLOC = AC2(IP,:,:)
               CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE.) 
               IMATDAA(IP,:,:) = IMATDA
               IMATRAA(IP,:,:) = IMATRA
             END IF !
           ELSE
             IF (LSOUBOUND) THEN ! Source terms on boundary ...
               IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
                 ACLOC = AC2(IP,:,:)
                 CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE.)
                 IMATDAA(IP,:,:) = IMATDA
                 IMATRAA(IP,:,:) = IMATRA
               ENDIF
             ENDIF
           ENDIF  
         END DO

#if defined ST41 || defined ST42
         LFIRSTSOURCE = .FALSE.
#endif
         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SOURCE_INT_IMP_WAM_PRE()
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER        :: IP, IS, ID
         REAL(rkind)    :: ACLOC(MSC,MDC), VEC2RAD
         REAL(rkind)    :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

!$OMP WORKSHARE
         IMATDAA = 0.
         IMATRAA = 0.
!$OMP END WORKSHARE

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(IP,ACLOC,IMATDA,IMATRA)
         DO IP = 1, MNP
           IF ((ABS(IOBP(IP)) .NE. 1 .AND. IOBP(IP) .NE. 3)) THEN
             IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
               THWOLD(:,1) = THWNEW
               U10NEW = MAX(TWO,SQRT(WINDXY(:,1)**2+WINDXY(:,2)**2)) * WINDFAC
               DO IS = 1, MSC
                 DO ID = 1, MDC
                   FL3(IP,ID,IS) = AC2(IP,IS,ID) * PI2 * SPSIG(IS)
                   FL(IP,ID,IS)  = FL3(IP,ID,IS)
                   SL(IP,ID,IS)  = FL(IP,ID,IS)
                 ENDDO
               END DO
               Z0NEW(IP) = Z0OLD(IP,1)
               THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))
               IF (LOUTWAM .AND. IP == TESTNODE) THEN
                 WRITE(111112,'(A10,I10)') 'BEFORE', IP
                 WRITE(111112,'(A10,F20.10)') 'FL3', SUM(FL3(IP,:,:))
                 WRITE(111112,'(A10,F20.10)') 'FL', SUM(FL(IP,:,:))
                 WRITE(111112,'(A10,F20.10)') 'THWOLD', THWOLD(IP,1)
                 WRITE(111112,'(A10,F20.10)') 'USOLD', USOLD(IP,1)
                 WRITE(111112,'(A10,F20.10)') 'U10NEW', U10NEW(IP)
                 WRITE(111112,'(A10,F20.10)') 'THWNEW', THWNEW(IP)
                 WRITE(111112,'(A10,F20.10)') 'Z0OLD', Z0OLD(IP,1)
                 WRITE(111112,'(A10,F20.10)') 'TAUW', TAUW(IP)
                 WRITE(111112,'(A10,F20.10)') 'ROAIRO', ROAIRO(IP,1)
                 WRITE(111112,'(A10,F20.10)') 'ZIDLOLD', ZIDLOLD(IP,1)
                 WRITE(111112,'(A10,F20.10)') 'Z0NEW', Z0NEW(IP)
                 WRITE(111112,'(A10,F20.10)') 'ROAIRN', ROAIRN(IP)
                 WRITE(111112,'(A10,F20.10)') 'ZIDLNEW', ZIDLNEW(IP)
                 WRITE(111112,'(A10,F20.10)') 'SL', SUM(SL(IP,:,:))
                 WRITE(111112,'(A10,F20.10)') 'FCONST', SUM(FCONST(IP,:))
               ENDIF
               CALL PREINTRHS (FL3(IP,:,:), FL(IP,:,:), IP, IP, 1, &
     &                         THWOLD(IP,1), USOLD(IP,1), &
     &                         TAUW(IP), Z0OLD(IP,1), &
     &                         ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                         U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                         Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                         SL(IP,:,:), FCONST(IP,:), FMEANWS(IP), MIJ(IP))
               IF (LOUTWAM .AND. IP == TESTNODE) THEN
                 WRITE(111112,'(A10,I10)') 'AFTER', IP
                 WRITE(111112,'(A10,F20.10)') 'FL3', SUM(FL3(IP,:,:))
                 WRITE(111112,'(A10,F20.10)') 'FL', SUM(FL(IP,:,:))
                 WRITE(111112,'(A10,F20.10)') 'THWOLD', THWOLD(IP,1)
                 WRITE(111112,'(A10,F20.10)') 'USOLD', USOLD(IP,1)
                 WRITE(111112,'(A10,F20.10)') 'U10NEW', U10NEW(IP)
                 WRITE(111112,'(A10,F20.10)') 'THWNEW', THWNEW(IP)
                 WRITE(111112,'(A10,F20.10)') 'Z0OLD', Z0OLD(IP,1)
                 WRITE(111112,'(A10,F20.10)') 'TAUW', TAUW(IP)
                 WRITE(111112,'(A10,F20.10)') 'ROAIRO', ROAIRO(IP,1)
                 WRITE(111112,'(A10,F20.10)') 'ZIDLOLD', ZIDLOLD(IP,1)
                 WRITE(111112,'(A10,F20.10)') 'Z0NEW', Z0NEW(IP)
                 WRITE(111112,'(A10,F20.10)') 'ROAIRN', ROAIRN(IP)
                 WRITE(111112,'(A10,F20.10)') 'ZIDLNEW', ZIDLNEW(IP)
                 WRITE(111112,'(A10,F20.10)') 'SL', SUM(SL(IP,:,:))
                 WRITE(111112,'(A10,F20.10)') 'FCONST', SUM(FCONST(IP,:))
               ENDIF
               DO IS = 1, MSC
                 DO ID = 1, MDC
                   IMATDAA(IP,IS,ID) =  FL(IP,ID,IS)
                   IMATRAA(IP,IS,ID) =  SL(IP,ID,IS)/PI2/SPSIG(IS)
                 ENDDO
               ENDDO 
               !ISELECT = 30
               !CALL SOURCETERMS(IP, AC2(IP,:,:), IMATRAA(IP,:,:), IMATDAA(IP,:,:), .FALSE.)
             END IF ! ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2)
           ELSE
             IF (LSOUBOUND) THEN ! Source terms on boundary ...
               IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
                 THWOLD(:,1) = THWNEW
                 U10NEW = MAX(TWO,SQRT(WINDXY(:,1)**2+WINDXY(:,2)**2)) * WINDFAC
                 DO IS = 1, MSC
                   DO ID = 1, MDC
                     FL3(IP,ID,IS) =  AC2(IP,IS,ID) * PI2 * SPSIG(IS)
                     FL(IP,ID,IS)  =  FL3(IP,ID,IS)
                     SL(IP,ID,IS)  =  FL(IP,ID,IS)
                   END DO
                 END DO
                 Z0NEW(IP) = Z0OLD(IP,1)
                 THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))
                 IF (LOUTWAM .AND. IP == TESTNODE) THEN
                   WRITE(111112,'(A10,I10)') 'BEFORE', IP
                   WRITE(111112,'(A10,F20.10)') 'FL3', SUM(FL3(IP,:,:))
                   WRITE(111112,'(A10,F20.10)') 'FL', SUM(FL(IP,:,:))
                   WRITE(111112,'(A10,F20.10)') 'THWOLD', THWOLD(IP,1)
                   WRITE(111112,'(A10,F20.10)') 'USOLD', USOLD(IP,1)
                   WRITE(111112,'(A10,F20.10)') 'U10NEW', U10NEW(IP)
                   WRITE(111112,'(A10,F20.10)') 'THWNEW', THWNEW(IP)
                   WRITE(111112,'(A10,F20.10)') 'Z0OLD', Z0OLD(IP,1)
                   WRITE(111112,'(A10,F20.10)') 'TAUW', TAUW(IP)
                   WRITE(111112,'(A10,F20.10)') 'ROAIRO', ROAIRO(IP,1)
                   WRITE(111112,'(A10,F20.10)') 'ZIDLOLD', ZIDLOLD(IP,1)
                   WRITE(111112,'(A10,F20.10)') 'Z0NEW', Z0NEW(IP)
                   WRITE(111112,'(A10,F20.10)') 'ROAIRN', ROAIRN(IP)
                   WRITE(111112,'(A10,F20.10)') 'ZIDLNEW', ZIDLNEW(IP)
                   WRITE(111112,'(A10,F20.10)') 'SL', SUM(SL(IP,:,:))
                   WRITE(111112,'(A10,F20.10)') 'FCONST', SUM(FCONST(1,:))
                 ENDIF
                 CALL PREINTRHS (FL3(IP,:,:), FL(IP,:,:), IP, IP, 1, &
     &                           THWOLD(IP,1), USOLD(IP,1), &
     &                           TAUW(IP), Z0OLD(IP,1), &
     &                           ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                           U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                           Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                           SL(IP,:,:), FCONST(IP,:), FMEANWS(IP), MIJ(IP))
                 IF (LOUTWAM .AND. IP == TESTNODE) THEN
                   WRITE(111112,'(A10,I10)') 'AFTER', IP
                   WRITE(111112,'(A10,F20.10)') 'FL3', SUM(FL3(IP,:,:))
                   WRITE(111112,'(A10,F20.10)') 'FL', SUM(FL(IP,:,:))
                   WRITE(111112,'(A10,F20.10)') 'THWOLD', THWOLD(IP,1)
                   WRITE(111112,'(A10,F20.10)') 'USOLD', USOLD(IP,1)
                   WRITE(111112,'(A10,F20.10)') 'U10NEW', U10NEW(IP)
                   WRITE(111112,'(A10,F20.10)') 'THWNEW', THWNEW(IP)
                   WRITE(111112,'(A10,F20.10)') 'Z0OLD', Z0OLD(IP,1)
                   WRITE(111112,'(A10,F20.10)') 'TAUW', TAUW(IP)
                   WRITE(111112,'(A10,F20.10)') 'ROAIRO', ROAIRO(IP,1)
                   WRITE(111112,'(A10,F20.10)') 'ZIDLOLD', ZIDLOLD(IP,1)
                   WRITE(111112,'(A10,F20.10)') 'Z0NEW', Z0NEW(IP)
                   WRITE(111112,'(A10,F20.10)') 'ROAIRN', ROAIRN(IP)
                   WRITE(111112,'(A10,F20.10)') 'ZIDLNEW', ZIDLNEW(IP)
                   WRITE(111112,'(A10,F20.10)') 'SL', SUM(SL(IP,:,:))
                   WRITE(111112,'(A10,F20.10)') 'FCONST', SUM(FCONST(IP,:))
                 ENDIF
                 DO IS = 1, MSC
                   DO ID = 1, MDC
                     IMATDAA(IP,IS,ID) = FL(IP,ID,IS)
                     IMATRAA(IP,IS,ID) = SL(IP,ID,IS)/PI2/SPSIG(IS)
                   ENDDO
                 ENDDO
                 !ISELECT = 30
                 !CALL SOURCETERMS(IP, AC2(IP,:,:), IMATRAA(IP,:,:), IMATDAA(IP,:,:), .FALSE.)
               ENDIF
             ENDIF
           ENDIF
         ENDDO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SOURCE_INT_IMP_WAM_POST
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER        :: IP, IS, ID
         REAL(rkind)    :: ACLOC(MSC,MDC), VEC2RAD
         REAL(rkind)    :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

!$OMP WORKSHARE
         IMATDAA = 0.
         IMATRAA = 0.
!$OMP END WORKSHARE

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(IP,ACLOC,IMATDA,IMATRA)
         DO IP = 1, MNP
           IF ((ABS(IOBP(IP)) .NE. 1 .AND. IOBP(IP) .NE. 3)) THEN
             IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
               THWOLD(:,1) = THWNEW
               U10NEW = MAX(TWO,SQRT(WINDXY(:,1)**2+WINDXY(:,2)**2)) * WINDFAC
               DO IS = 1, MSC
                 DO ID = 1, MDC
                   FL3(IP,ID,IS) =  AC2(IP,IS,ID) * PI2 * SPSIG(IS)
                   FL(IP,ID,IS)  =  IMATDAA(IP,IS,ID)
                   SL(IP,ID,IS)  =  IMATRAA(IP,IS,ID) * PI2 * SPSIG(IS)
                 END DO
                 Z0NEW(IP) = Z0OLD(IP,1)
               END DO
               THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))
               IF (IP == TESTNODE) WRITE(*,'(A20,3F15.8)') 'POST BEFORE', SUM(SL(IP,:,:)), SUM(FL3(IP,:,:)),  SUM(FL(IP,:,:))
               CALL POSTINTRHS (FL3(IP,:,:), FL(IP,:,:), IP, IP, 1, &
     &                          THWOLD(IP,1), USOLD(IP,1), &
     &                          TAUW(IP), Z0OLD(IP,1), &
     &                          ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                          U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                          Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                          SL(IP,:,:), FCONST(IP,:), FMEANWS(IP), MIJ(IP))
               IF (IP == TESTNODE) WRITE(*,'(A20,3F15.8)') 'POST AFTER', SUM(SL(IP,:,:)), SUM(FL3(IP,:,:)),  SUM(FL(IP,:,:))
               DO IS = 1, MSC
                 DO ID = 1, MDC
                   IMATDAA(IP,IS,ID) = FL(IP,ID,IS)
                   IMATRAA(IP,IS,ID) = SL(IP,ID,IS)/PI2/SPSIG(IS)
                 ENDDO
               ENDDO
             END IF !
           ELSE
             IF (LSOUBOUND) THEN ! Source terms on boundary ...
               IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
                 FL = FL3
                 THWOLD(:,1) = THWNEW
                 U10NEW = MAX(TWO,SQRT(WINDXY(:,1)**2+WINDXY(:,2)**2)) * WINDFAC
                 DO IS = 1, MSC
                   DO ID = 1, MDC
                     FL3(IP,ID,IS) =  AC2(IP,IS,ID) * PI2 * SPSIG(IS)
                     FL(IP,ID,IS)  =  IMATDAA(IP,IS,ID)
                     SL(IP,ID,IS)  =  IMATRAA(IP,IS,ID) * PI2 * SPSIG(IS)
                   END DO
                   Z0NEW(IP) = Z0OLD(IP,1)
                 END DO
                 THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))
                 CALL POSTINTRHS (FL3(IP,:,:), FL(IP,:,:), IP, IP, 1, &
     &                            THWOLD(IP,1), USOLD(IP,1), &
     &                            TAUW(IP), Z0OLD(IP,1), &
     &                            ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                            U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                            Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                            SL(IP,:,:), FCONST(IP,:), FMEANWS(IP), MIJ(IP))
                 DO IS = 1, MSC
                   DO ID = 1, MDC
                     IMATDAA(IP,IS,ID) = FL(IP,ID,IS)
                     IMATRAA(IP,IS,ID) = SL(IP,ID,IS)/PI2/SPSIG(IS)
                   ENDDO
                 ENDDO
               ENDIF
             ENDIF
           ENDIF
           DO IS = 1, MSC
             DO ID = 1, MDC
               AC2(IP,IS,ID) = FL3(IP,ID,IS) / PI2 / SPSIG(IS)
             END DO
           END DO
         END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INT_IP_STAT(IP,DT,LIMITER,ACLOC)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP
         LOGICAL, INTENT(IN) :: LIMITER
         REAL(rkind), INTENT(IN) :: DT
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)
         REAL(rkind)                :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER :: IS, ID

         REAL(rkind)    :: NPF(MSC)
         REAL(rkind)    :: OLDAC
         REAL(rkind)    :: NEWDAC
         REAL(rkind)    :: MAXDAC, CONST, SND, USTAR

         ISELECT = 10
         CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE.)  ! 1. CALL

         CONST = PI2**2*3.0*1.0E-7*DT*SPSIG(MSC_HF(IP))
         SND   = PI2*5.6*1.0E-3

         IF (LIMITER) THEN
           IF (MELIM == 1) THEN
             NPF = 0.0081_rkind*LIMFAK/(TWO*SPSIG*WK(IP,:)**3*CG(IP,:))
           ELSE IF (MELIM == 2) THEN
             DO IS = 1, MSC
               USTAR = MAX(UFRIC(IP), G9*SND/SPSIG(IS))
               NPF(IS) = ABS((CONST*USTAR)/(SPSIG(IS)**3*WK(IP,IS)))
             END DO
           END IF
         END IF

!         if (SUM(ACLOC) .NE. SUM(ACLOC)) STOP 'NAN l. 174 wwm_specint.F90'

         DO IS = 1, MSC
           MAXDAC = NPF(IS)
           DO ID = 1, MDC
             OLDAC  = ACLOC(IS,ID)
             NEWDAC = IMATRA(IS,ID) * DT / ( 1.0 - DT * MIN(ZERO,IMATDA(IS,ID))) 
             IF (LIMITER) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
             ACLOC(IS,ID) = MAX( 0.0_rkind, OLDAC + NEWDAC ) 
           END DO
         END DO

!         write(*,'(A10,2I10,L10,I10,2F15.6)') 'after', ip, iobp(ip), limiter, iselect, sum(acloc), sum(imatra)

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE RKS_SP3(IP,DTSII,LIMITER,ACLOC)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP
         LOGICAL, INTENT(IN) :: LIMITER
         REAL(rkind), INTENT(IN)    :: DTSII
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)

         INTEGER :: IS, ID
         REAL(rkind)    :: ACOLD(MSC,MDC)
         REAL(rkind)    :: NPF(MSC)
         REAL(rkind)    :: NEWDAC, MAXDAC, CONST, SND, USTAR
         REAL(rkind)    :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         CONST = PI2**2*3.0*1.0E-7*DTSII*SPSIG(MSC_HF(IP))
         SND   = PI2*5.6*1.0E-3

         IF (LIMITER) THEN
           IF (MELIM == 1) THEN
             NPF = 0.0081_rkind*LIMFAK/(TWO*SPSIG*WK(IP,:)**3*CG(IP,:))
           ELSE IF (MELIM == 2) THEN
             DO IS = 1, MSC
               USTAR = MAX(UFRIC(IP), G9*SND/SPSIG(IS))
               NPF(IS) = ABS((CONST*USTAR)/(SPSIG(IS)**3*WK(IP,IS)))
             END DO
           END IF
         END IF

         CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE.)  ! 1. CALL

         ACOLD = ACLOC

         DO IS = 1, MSC
           MAXDAC = NPF(IS)
           DO ID = 1, MDC
             NEWDAC = IMATRA(IS,ID) * DTSII / ( ONE - DTSII * MIN(ZERO,IMATDA(IS,ID)) )
             IF (LIMITER) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
             ACLOC(IS,ID) = MAX( ZERO, ACLOC(IS,ID) + NEWDAC )
           END DO
         END DO

         !WRITE(*,*) '1 RK-TVD', SUM(ACOLD), SUM(ACLOC)

         CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE.) ! 2. CALL

         DO IS = 1, MSC
           MAXDAC = NPF(IS)
           DO ID = 1, MDC
             NEWDAC = IMATRA(IS,ID) * DTSII / ( ONE - DTSII * MIN(ZERO,IMATDA(IS,ID)) )
             IF (LIMITER) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC) 
             ACLOC(IS,ID) = MAX( ZERO, 0.75_rkind * ACOLD(IS,ID) +  0.25_rkind * ACLOC(IS,ID) + 0.25_rkind * NEWDAC)
           END DO
         END DO

         !WRITE(*,*) '2 RK-TVD', SUM(ACOLD), SUM(ACLOC)

         CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE.) ! 3. CALL

         DO IS = 1, MSC
           MAXDAC = NPF(IS)
           DO ID = 1, MDC
             NEWDAC = IMATRA(IS,ID) * DTSII / ( ONE - DTSII * MIN(ZERO,IMATDA(IS,ID)) )
             IF (LIMITER) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
             ACLOC(IS,ID) = MAX( ZERO, ONE/THREE * ACOLD(IS,ID) + TWO/THREE * ACLOC(IS,ID) + TWO/THREE * NEWDAC)
           END DO
         END DO

         !WRITE(*,*) '3 RK-TVD', SUM(ACOLD), SUM(ACLOC)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INT_IP_DYN(IP, DT, LIMIT, DTMIN, ITRMX, ACLOC, ITER)

         USE DATAPOOL
#ifdef ST41
         USE W3SRC4MD_OLD
#elif ST42
         USE W3SRC4MD
#endif
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IP, ITRMX
         INTEGER, INTENT(OUT)       :: ITER
         LOGICAL,INTENT(IN)         :: LIMIT
         REAL(rkind), INTENT(IN)    :: DTMIN, DT
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)
         REAL(rkind)                :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER :: IS, ID

         REAL(rkind)    :: ACOLD(MSC,MDC)
         REAL(rkind)    :: NPF(MSC)
         REAL(rkind)    :: AFILT
         REAL(rkind)    :: NEWDAC, CONST, SND, DAMAX, AFAC
         REAL(rkind)    :: MAXDAC, DTMAX, DTTOT, DTLEFT, DT4SI, USTAR

         REAL(rkind), PARAMETER :: MAXDTFAC = VERYLARGE 

         LOGICAL :: LSTABLE

         REAL(rkind) :: XPP, XRR, XFILT, XREL, XFLT, FACP, DAM(MSC)

#ifdef ST_DEF
         XPP    = 0.15
         XRR    = 0.10
         XFILT  = 0.0001 ! AR: check why it must be so small ..
         XPP    = MAX ( 1.E-6_rkind , XPP )
         XRR    = MAX ( 1.E-6_rkind , XRR )
         XREL   = XRR
         XFLT   = MAX ( ZERO , XFILT ) 
         FACP   = 2*XPP / PI2 * 0.62E-3 * PI2**4 / G9**2  ! s4/m4
         DAM    = FACP / ( SPSIG * WK(IP,:)**3 ) / CG(IP,:) ! s * m³ * s4/m4 = 
         AFILT  = MAX ( DAM(MSC) , XFLT*MAXVAL(ACLOC))!*PI2*SPSIG(MSC_HF(IP)) )
#endif

         CONST = PI2**2*3.0*1.0E-7*DT*SPSIG(MSC_HF(IP))
         SND   = PI2*5.6*1.0E-3

         IF (LIMIT) THEN
           IF (MELIM == 1) THEN
             NPF = 0.0081_rkind*LIMFAK/(TWO*SPSIG*WK(IP,:)**3*CG(IP,:))
           ELSE IF (MELIM == 2) THEN
             DO IS = 1, MSC
               USTAR = MAX(UFRIC(IP), G9*SND/SPSIG(IS))
               NPF(IS) = ABS((CONST*USTAR)/(SPSIG(IS)**3*WK(IP,IS)))
             END DO
           END IF
         END IF

         DTTOT = 0.
         ITER  = 0

         DO WHILE ( DTTOT < DT )

           ACOLD = ACLOC
           IF (ITER == 0) DT4SI = DT
           IF (ITER == 0) DTLEFT = DT
           ITER  = ITER + 1
           DTMAX =  DT

           CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE.)  ! 1. CALL
        
           DO ID = 1, MDC
             DO IS = 1, MSC_HF(IP)
               IF (ABS(IMATRA(IS,ID)) .GT. VERYSMALL) THEN                
                 DAMAX  = MIN ( DAM(IS) , MAX ( XREL*ACLOC(IS,ID), AFILT ) )
                 AFAC  = ONE / MAX( 1.E-10_rkind , ABS(IMATRA(IS,ID)/DAMAX) )
                 DTMAX = MIN (DTMAX ,AFAC/(MAX(1.E-10_rkind, ONE + AFAC*MIN(ZERO,IMATDA(IS,ID))))) 
               END IF
             END DO
           END DO

           DT4SI  = DTMAX
           DT4SI  = MAX(DTMIN,DT4SI) ! This is not entirely stable !!!
           DTLEFT = DT - DTTOT

           IF ( DTLEFT > THR .AND. DTLEFT < DT4SI) THEN
             DT4SI = (DT - DTTOT)
           ELSE IF ( DTLEFT .GE. DT4SI .AND. ITER .EQ. ITRMX) THEN
             DT4SI = DTLEFT 
             LSTABLE = .FALSE. 
           END IF 

           DTTOT = DTTOT + DT4SI

           IF ( DT4SI .LT. DTMAX + SMALL) THEN
             LSTABLE = .TRUE. 
           ELSE
             LSTABLE = .FALSE.
           END IF

           IF (LSTABLE) THEN
             DO IS = 1, MSC_HF(IP)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(ZERO,IMATDA(IS,ID)) )
                 ACLOC(IS,ID) = MAX( ZERO, ACLOC(IS,ID) + NEWDAC )
               END DO
             END DO
             CALL SOURCETERMS(IP,ACLOC,IMATRA,IMATDA,.FALSE.) ! 2. CALL
             DO IS = 1, MSC_HF(IP)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(ZERO,IMATDA(IS,ID)) )
                 ACLOC(IS,ID) = MAX( ZERO, 3._rkind/4._rkind * ACOLD(IS,ID) + 1._rkind/4._rkind * ACLOC(IS,ID) + 1._rkind/4._rkind * NEWDAC)
               END DO
             END DO
             CALL SOURCETERMS(IP,ACLOC,IMATRA,IMATDA, .FALSE.) ! 3. CALL
             DO IS = 1, MSC_HF(IP)
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(ZERO,IMATDA(IS,ID)) )
                 ACLOC(IS,ID) = MAX( ZERO,  1._rkind/3._rkind * ACOLD(IS,ID) + 2._rkind/3._rkind * ACLOC(IS,ID) + 2._rkind/3._rkind * NEWDAC)
               END DO
             END DO
           ELSE ! .NOT. LSTABLE
             DO IS = 1, MSC_HF(IP)
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(0._rkind,IMATDA(IS,ID)) )
                 IF (LIMIT) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
                 ACLOC(IS,ID) = MAX( ZERO, ACLOC(IS,ID) + NEWDAC )
               END DO
             END DO
             CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE.) ! 2. CALL
             DO IS = 1, MSC_HF(IP)
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(0._rkind,IMATDA(IS,ID)) )
                 IF (LIMIT) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
                 ACLOC(IS,ID) = MAX( ZERO, 3./4. * ACOLD(IS,ID) +  1./4. * ACLOC(IS,ID) + 1./4. * NEWDAC)
               END DO
             END DO
             CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE.) ! 3. CALL
             DO IS = 1, MSC_HF(IP)
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(0._rkind,IMATDA(IS,ID)) )
                 IF (LIMIT) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
                 ACLOC(IS,ID) = MAX( ZERO, 1./3. * ACOLD(IS,ID) +  2./3. * ACLOC(IS,ID) + 2./3. * NEWDAC)
               END DO
             END DO
           END IF
         END DO      

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
