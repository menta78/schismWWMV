#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
#define ACTIVATE_SMETHOD_5
#undef ACTIVATE_SMETHOD_5
      SUBROUTINE SOURCE_INT_EXP

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

!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP&         SHARED(MNP,MSC,MDC,DEP,DMIN,IOBP,SMETHOD, &
!$OMP&         DT4S,DTMIN_SDS,DTMIN_SIN,DTMIN_SBR,DTMIN_DYN, &
!$OMP&         DTMIN_SNL3, DTMIN_SNL4, DTMIN_SBF, &
!$OMP&         LSOUBOUND,ISHALLOW,LADVTEST,LMAXETOT, &
!$OMP&         NDYNITER_SIN,NDYNITER_SNL4, NDYNITER_SDS, &
!$OMP&         NDYNITER_SBR, NDYNITER_SNL3, NDYNITER_SBF, &
!$OMP&         NDYNITER,LSOURCESWAM,LLIMT) &
!$OMP&         PRIVATE(IP,IS,ID,ACLOC,AC2,AC1, NIT_SIN,NIT_SDS,&
!$OMP&         NIT_SNL4,NIT_SNL3,NIT_SBR,NIT_SBF,NIT_ALL,ISELECT,SSBRL2)
!$OMP DO SCHEDULE(DYNAMIC)





         DO IP = 1, MNP
!           IF (IP_IS_STEADY(IP) .EQ. 1) CYCLE
           IF (DEP(IP) .LT. DMIN) CYCLE
           IF (IOBP(IP) .EQ. 0) THEN
               ACLOC  = AC2(IP,:,:)
               IF (SMETHOD == 1) THEN
                 ISELECT = 30
                 CALL RKS_SP3(IP,DT4S,.FALSE.,ACLOC)
                 IF (LSOURCESWAM) THEN
                   CALL SOURCE_INT_EXP_WAM(IP, ACLOC)  
                 ELSE
                   ISELECT = 20
                   CALL INT_IP_STAT(IP,DT4S,LLIMT,ACLOC)
                 ENDIF
                 ISELECT = 4
                 CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_DYN, NDYNITER, ACLOC, NIT_ALL)
               ELSE IF (SMETHOD == 2) THEN
                 ISELECT = 10
                 CALL INT_IP_STAT(IP,DT4S, LLIMT,ACLOC)
               ELSE IF (SMETHOD == 3) THEN
                 ISELECT = 10
                 CALL RKS_SP3(IP,DT4S,LLIMT,ACLOC)
               ELSE IF (SMETHOD == 4) THEN
                 ISELECT = 10
                 CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_DYN, NDYNITER, ACLOC, NIT_ALL)
               ELSE IF (SMETHOD == 5) THEN ! Full splitting of all source embedded within a dynamic RK-3 Integration ... 
                 ISELECT = 1
                 CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SIN,  NDYNITER_SIN  , ACLOC, NIT_SIN) ! Sin
                 ISELECT = 2 
                 CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SNL4,  NDYNITER_SNL4 , ACLOC, NIT_SNL4)! Snl4b
                 ISELECT = 3 
                 CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SDS,  NDYNITER_SDS  , ACLOC, NIT_SDS) ! Sds
                 ISELECT = 4 
                 CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SNL3,  NDYNITER_SNL3 , ACLOC, NIT_SNL3)! Snl3
                 ISELECT = 5 
                 CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SBR,  NDYNITER_SBR  , ACLOC, NIT_SBR) ! Sbr
                 ISELECT = 6 
                 CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SBF,  NDYNITER_SBF  , ACLOC, NIT_SBF) ! Sbf
               END IF
               !ISELECT = 1
               !CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .TRUE.) ! Update everything based on the new spectrum ...
               IF (LMAXETOT .AND. .NOT. LADVTEST .AND. ISHALLOW(IP) .EQ. 1) THEN
                 CALL BREAK_LIMIT(IP,ACLOC,SSBRL2)
               ENDIF
               AC2(IP,:,:) = ACLOC
           ELSE !Boundary node ... 
             IF (LSOUBOUND) THEN ! Source terms on boundary ...
               IF (IOBP(IP) .NE. 2) THEN
                 ACLOC  = AC2(IP,:,:)
                 IF (SMETHOD == 1) THEN
                   ISELECT = 30
                   CALL RKS_SP3(IP,DT4S,.FALSE.,ACLOC)
                   IF (LSOURCESWAM) THEN
                     CALL SOURCE_INT_EXP_WAM(IP, ACLOC)
                   ELSE
                     ISELECT = 20
                     CALL INT_IP_STAT(IP,DT4S,LLIMT,ACLOC)
                   ENDIF
                   ISELECT = 4
                   CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_DYN, NDYNITER, ACLOC, NIT_ALL)
                 ELSE IF (SMETHOD == 2) THEN
                   ISELECT = 10
                   CALL INT_IP_STAT(IP,DT4S, LLIMT,ACLOC)
                 ELSE IF (SMETHOD == 3) THEN
                   ISELECT = 10
                   CALL RKS_SP3(IP,DT4S,LLIMT,ACLOC)
                 ELSE IF (SMETHOD == 4) THEN
                   ISELECT = 10
                   CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_DYN, NDYNITER, ACLOC, NIT_ALL)
                 ELSE IF (SMETHOD == 5) THEN ! Full splitting of all source embedded within a dynamic RK-3 Integration ... 
                   ISELECT = 1
                   CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SIN,  NDYNITER_SIN  , ACLOC, NIT_SIN) ! Sin
                   ISELECT = 2
                   CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SNL4,  NDYNITER_SNL4 , ACLOC, NIT_SNL4)! Snl4b
                   ISELECT = 3
                   CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SDS,  NDYNITER_SDS  , ACLOC, NIT_SDS) ! Sds
                   ISELECT = 4
                   CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SNL3,  NDYNITER_SNL3 , ACLOC, NIT_SNL3)! Snl3
                   ISELECT = 5
                   CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SBR,  NDYNITER_SBR  , ACLOC, NIT_SBR) ! Sbr
                   ISELECT = 6
                   CALL INT_IP_DYN(IP, DT4S, LLIMT, DTMIN_SBF,  NDYNITER_SBF  , ACLOC, NIT_SBF) ! Sbf
                 END IF
                 !ISELECT = 1
                 !CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .TRUE.) ! Update everything based on the new spectrum ...
                 IF (LMAXETOT .AND. .NOT. LADVTEST .AND. ISHALLOW(IP) .EQ. 1) THEN
                   CALL BREAK_LIMIT(IP,ACLOC,SSBRL2)
                 ENDIF
                 AC2(IP,:,:) = ACLOC
               ENDIF
             ELSE
               ACLOC = AC2(IP,:,:)
               IF (LMAXETOT .AND. .NOT. LADVTEST .AND. ISHALLOW(IP) .EQ. 1) THEN ! limit wave height on the boundary ...
                 CALL BREAK_LIMIT(IP,ACLOC,SSBRL2)
                 AC2(IP,:,:) = ACLOC
               ENDIF
             ENDIF
           ENDIF
           AC1(IP,:,:) = AC2(IP,:,:)
         ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL 

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
               !CALL CYCLE3(IP, ACLOC, IMATRA, IMATDA)
               IMATDAA(IP,:,:) = IMATDA 
               IMATRAA(IP,:,:) = IMATRA 
             ENDIF
           ELSE
             IF (LSOUBOUND) THEN ! Source terms on boundary ...
               IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
                 ACLOC = AC2(IP,:,:)
                 CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE.)
                 !CALL CYCLE3(IP, ACLOC, IMATRA, IMATDA)
                 IMATDAA(IP,:,:) = IMATDA
                 IMATRAA(IP,:,:) = IMATRA 
               ENDIF
             ENDIF
           ENDIF
         END DO
#if defined ST41 || defined ST42
         LFIRSTSOURCE = .FALSE.
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************      
      SUBROUTINE SOURCE_INT_EXP_WAM(IP, ACLOC)
         USE DATAPOOL
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: IP
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)
         
         INTEGER      :: IS, ID, IMETHOD
         REAL(rkind)  :: VEC2RAD
         REAL(rkind)  :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind)  :: ETOT,SME01,SME10,KME01,KMWAM,KMWAM2,HS,WIND10
         REAL(rkind)  :: ETAIL,EFTAIL,EMAX,LIMAC,NEWDAC,FPM,WINDTH,TEMP,GTEMP1
         REAL(rkind)  :: RATIO,LIMFAC,LIMDAC,GTEMP2,FLHAB,DELFL,USFM, NEWDACDT
         REAL(rkind)  :: MAXDAC, MAXDACDT, MAXDACDTDA, SC, SP, DNEWDACDTDA, JAC, FF

         REAL(rkind),DIMENSION(MDC,MSC)  :: SSDS,DSSDS,SSNL4,DSSNL4,SSIN,DSSIN
          
         IF (MESIN .EQ. 0 .AND. MESDS .AND. 0 .AND. MESNL .EQ. 0) THEN
           DO IS = 1, MSC
             DO ID = 1, MDC
               FL3(IP,ID,IS) = ACLOC(IS,ID) * PI2 * SPSIG(IS)
               FL(IP,ID,IS)  = FL3(IP,ID,IS)
               SL(IP,ID,IS)  = FL(IP,ID,IS)
             END DO
           END DO
           THWOLD(IP,1) = THWNEW(IP)
           U10NEW(IP) = MAX(TWO,SQRT(WINDXY(IP,1)**2+WINDXY(IP,2)**2))*WINDFAC
           Z0NEW(IP) = Z0OLD(IP,1)
           THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))
           IF (.FALSE.) THEN
             CALL IMPLSCH (FL3(1,:,:), FL(1,:,:), IP, IP, 1, &
     &                     THWOLD(IP,1), USOLD(IP,1), &
     &                     TAUW(IP), Z0OLD(IP,1), &
     &                     ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                     U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                     Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                     SL(1,:,:), FCONST(1,:))
           ELSE
             CALL PREINTRHS (FL3(1,:,:), FL(1,:,:), IP, IP, 1, &
     &                       THWOLD(IP,1), USOLD(IP,1), &
     &                       TAUW(IP), Z0OLD(IP,1), &
     &                       ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                       U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                       Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                       SL(1,:,:), FCONST(1,:), FMEANWS(IP), MIJ(IP), &
     &                       SSDS, DSSDS, SSIN, DSSIN, &
     &                       SSNL4, DSSNL4)
             CALL INTSPECWAM (FL3(1,:,:), FL(1,:,:), IP, IP, 1, &
     &                      THWOLD(IP,1), USOLD(IP,1), &
     &                      TAUW(IP), Z0OLD(IP,1), &
     &                      ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                      U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                      Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                      SL(1,:,:), FCONST(1,:), FMEANWS(IP), MIJ(IP))
             CALL POSTINTRHS (FL3(1,:,:), FL(1,:,:), IP, IP, 1, &
     &                      THWOLD(IP,1), USOLD(IP,1), &
     &                      TAUW(IP), Z0OLD(IP,1), &
     &                      ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                      U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                      Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                      SL(1,:,:), FCONST(1,:), FMEANWS(IP), MIJ(IP))
           ENDIF ! true false ...
           DO IS = 1, MSC
             DO ID = 1, MDC
               ACLOC(IS,ID) =  FL3(1,ID,IS) / PI2 / SPSIG(IS)
             END DO
           END DO
         ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SOURCE_INT_IMP_WAM_PRE
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER      :: IP, IS, ID, IMETHOD
         REAL(rkind)  :: ACLOC(MSC,MDC), VEC2RAD
         REAL(rkind)  :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind)  :: SSNL3(MSC,MDC),DSSNL3(MSC,MDC)
         REAL(rkind)  :: SSBR(MSC,MDC),DSSBR(MSC,MDC)
         REAL(rkind)  :: SSBF(MSC,MDC),DSSBF(MSC,MDC), SSINL(MSC,MDC)
         REAL(rkind)  :: ETOT,SME01,SME10,KME01,KMWAM,KMWAM2,HS,WIND10
         REAL(rkind)  :: ETAIL,EFTAIL,EMAX,LIMAC,NEWDAC,FPM,WINDTH,TEMP,GTEMP1
         REAL(rkind)  :: RATIO,LIMFAC,LIMDAC,GTEMP2,FLHAB,DELFL,USFM, NEWDACDT
         REAL(rkind)  :: MAXDAC, MAXDACDT, MAXDACDTDA, SC, SP, DNEWDACDTDA, JAC, FF

         REAL(rkind),DIMENSION(MDC,MSC)  :: SSDS,DSSDS,SSNL4,DSSNL4,SSIN,DSSIN

         IMETHOD = 4 

!$OMP WORKSHARE
         IMATDAA = ZERO
         IMATRAA = ZERO
!$OMP END WORKSHARE

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(IP,ACLOC,IMATDA,IMATRA)
         DO IP = 1, MNP
           SSNL3 = ZERO; DSSNL3 = ZERO
           SSBR  = ZERO; DSSBR  = ZERO
           SSBF  = ZERO; DSSBF  = ZERO
           SSINL = ZERO
           IF (DEP(IP) .LT. DMIN) THEN
             IMATRAA(IP,:,:) = ZERO
             IMATDAA(IP,:,:) = ZERO
             CYCLE
           ENDIF
           IF (MESIN .EQ. 0 .AND. MESDS .AND. 0 .AND. MESNL .EQ. 0) THEN
             DO IS = 1, MSC
               DO ID = 1, MDC
                 FL3(IP,ID,IS) = AC2(IP,IS,ID) * PI2 * SPSIG(IS)
                 FL(IP,ID,IS)  = FL3(IP,ID,IS)
                 SL(IP,ID,IS)  = FL(IP,ID,IS)
               END DO
             END DO
             THWOLD(IP,1) = THWNEW(IP)
             U10NEW(IP) = MAX(TWO,SQRT(WINDXY(IP,1)**2+WINDXY(IP,2)**2))*WINDFAC
             Z0NEW(IP) = Z0OLD(IP,1)
             THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))
             CALL PREINTRHS (FL3(IP,:,:), FL(IP,:,:), IP, IP, 1, &
     &                       THWOLD(IP,1), USOLD(IP,1), &
     &                       TAUW(IP), Z0OLD(IP,1), &
     &                       ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                       U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                       Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                       SL(IP,:,:), FCONST(IP,:), FMEANWS(IP), MIJ(IP), &
     &                       SSDS, DSSDS, SSIN, DSSIN, &
     &                       SSNL4, DSSNL4)
           ENDIF ! MESIN .EQ. 0 .AND. MESDS .AND. 0 .AND. MESNL = 0
           IF (IOBP(IP) .EQ. 0) THEN
             DO ID = 1, MDC
               DO IS = 1, MSC 
                 IF (AC2(IP,IS,ID) .LT. THR) CYCLE
                 JAC = ONE/PI2/SPSIG(IS)
                 !IMATDAA(IP,IS,ID) =  FL(IP,ID,IS) !... this is not working right, reason is unknown, there signchanges that should not be
                 !IMATRAA(IP,IS,ID) =  SL(IP,ID,IS)/PI2/SPSIG(IS) 
                 FF = FL3(IP,ID,IS)
                 !WRITE(11140,'(2I10,7E20.10)') IS, ID, FF, SSDS(ID,IS)/FF, DSSDS(ID,IS), SSIN(ID,IS)/FF, DSSNL4(ID,IS), SSNL4(ID,IS)/FF, DSSNL4(ID,IS) 
                 IF (IMETHOD == 0) THEN
                   IMATRAA(IP,IS,ID) =  SL(IP,ID,IS)/PI2/SPSIG(IS)
                   IMATDAA(IP,IS,ID) = ZERO 
                 ELSE IF (IMETHOD == 1) THEN 
                   IMATRAA(IP,IS,ID) = (SSIN(ID,IS)+SSDS(ID,IS)+SSNL4(ID,IS))*JAC
                 ELSE IF (IMETHOD == 2) THEN
                   IMATRAA(IP,IS,ID) = (SSIN(ID,IS)+SSNL4(ID,IS))*JAC
                   IMATDAA(IP,IS,ID) = -TWO*DSSDS(ID,IS)
                 ELSE IF (IMETHOD == 3) THEN
                   IF (SSIN(ID,IS) .GT. ZERO) THEN
                     IMATRAA(IP,IS,ID) = SSIN(ID,IS)*JAC
                   ELSE
                     IMATDAA(IP,IS,ID) = -TWO*DSSIN(ID,IS)
                   ENDIF
                 ELSE IF (IMETHOD == 4) THEN
                   GTEMP1 = MAX((1.-DT4A*FL(IP,ID,IS)),1.)
                   GTEMP2 = DT4A*SL(IP,ID,IS)/GTEMP1
                   FLHAB  = ABS(GTEMP2)
                   DELFL  = COFRM4(IS)*DT4A
                   USFM   = USNEW(IP)*MAX(FMEANWS(IP),FMEAN(IP))
                   TEMP   = USFM*DELFL
                   FLHAB  = MIN(FLHAB,TEMP)
                   IMATRAA(IP,IS,ID) = SIGN(FLHAB,GTEMP2)/DT4A*JAC
                   LIMFAC            = MIN(ONE,ABS(SIGN(FLHAB,GTEMP2/DT4A))/MAX(THR,ABS(IMATRAA(IP,IS,ID))))
                   IMATDAA(IP,IS,ID) = ZERO ! -FL(IP,ID,IS)
                 ENDIF 
               ENDDO ! ID
             ENDDO ! IS
             ACLOC = AC2(IP,:,:)
             IF (.NOT. LINID) THEN
               CALL SET_WIND( IP, WIND10, WINDTH )
               CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
               CALL SIN_LIN_CAV(IP,WINDTH,FPM,IMATRA,SSINL)
               IMATRAA(IP,:,:) = IMATRAA(IP,:,:) + SSINL
             ENDIF
             IF (ISHALLOW(IP) .EQ. 1) THEN
               CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)
               IF (MESTR .GT. 0) THEN
                 CALL TRIADSWAN_NEW (IP, HS, SME01, ACLOC, IMATRA, IMATDA, SSNL3, DSSNL3)
                 DO IS = 1, MSC
                   DO ID = 1, MDC
                     NEWDAC = SSNL3(IS,ID)*DT4A/MAX((1.-DT4A*DSSNL3(IS,ID)),1.)
                     MAXDAC = 0.0081*LIMFAK/(TWO*SPSIG(IS)*WK(IP,IS)**3*CG(IP,IS))*100
                     LIMFAC = ONE/MAX(ONE,NEWDAC/MAXDAC)
                     SC = SIGN(MIN(ABS(NEWDAC),MAXDAC),NEWDAC)/DT4A
                     !SSNL3(IS,ID)  = SC
                     !DSSNL3(IS,ID) = DSSNL3(IS,ID)*LIMFAC
                     !IF (ABS(SC) .GT. THR) WRITE(*,'(2I10,5F20.8)') IS, ID, NEWDAC, MAXDAC, DSSNL3(IS,ID), LIMFAC
                   END DO
                 END DO
               ENDIF ! MESTR
               IF (MESBR .GT. 0) CALL SDS_SWB(IP, SME01, KMWAM, ETOT, HS, ACLOC, IMATRA, IMATDA, SSBR, DSSBR)
               IF (MESBF .GT. 0) CALL SDS_BOTF(IP,ACLOC,IMATRA,IMATDA,SSBF,DSSBF)
               IMATDAA(IP,:,:) = IMATDAA(IP,:,:) + DSSBR  + DSSNL3 + DSSBF
               IMATRAA(IP,:,:) = IMATRAA(IP,:,:) + SSBR + SSNL3 
             ENDIF ! ISHALLOW(IP) .EQ. 1
           ELSE ! IOBP(IP) .NE. 0
             IF (LSOUBOUND) THEN ! Source terms on boundary ...
               IF (IOBP(IP) .NE. 2) THEN
                 DO ID = 1, MDC
                   DO IS = 1, MSC
                     IF (AC2(IP,IS,ID) .LT. THR) CYCLE
                     JAC = ONE/PI2/SPSIG(IS)
                     !IMATDAA(IP,IS,ID) =  FL(IP,ID,IS) !... this is not working right, reason is unknown, there signchanges that should not be
                     !IMATRAA(IP,IS,ID) =  SL(IP,ID,IS)/PI2/SPSIG(IS) 
                     FF = FL3(IP,ID,IS)
                     !WRITE(11140,'(2I10,7E20.10)') IS, ID, FF, SSDS(ID,IS)/FF, DSSDS(ID,IS), SSIN(ID,IS)/FF, DSSNL4(ID,IS), SSNL4(ID,IS)/FF, DSSNL4(ID,IS) 
                     IF (IMETHOD == 0) THEN
                       IMATRAA(IP,IS,ID) =  SL(IP,ID,IS)/PI2/SPSIG(IS)
                       IMATDAA(IP,IS,ID) = ZERO
                     ELSE IF (IMETHOD == 1) THEN
                       IMATRAA(IP,IS,ID) = (SSIN(ID,IS)+SSDS(ID,IS)+SSNL4(ID,IS))*JAC
                     ELSE IF (IMETHOD == 2) THEN
                       IMATRAA(IP,IS,ID) = (SSIN(ID,IS)+SSNL4(ID,IS))*JAC
                       IMATDAA(IP,IS,ID) = -TWO*DSSDS(ID,IS)
                     ELSE IF (IMETHOD == 3) THEN
                       IF (SSIN(ID,IS) .GT. ZERO) THEN
                         IMATRAA(IP,IS,ID) = SSIN(ID,IS)*JAC
                       ELSE
                         IMATDAA(IP,IS,ID) = -TWO*DSSIN(ID,IS)
                       ENDIF
                     ELSE IF (IMETHOD == 4) THEN
                       GTEMP1 = MAX((1.-DT4A*FL(IP,ID,IS)),1.)
                       GTEMP2 = DT4A*SL(IP,ID,IS)/GTEMP1
                       FLHAB  = ABS(GTEMP2)
                       DELFL  = COFRM4(IS)*DT4A
                       USFM   = USNEW(IP)*MAX(FMEANWS(IP),FMEAN(IP))
                       TEMP   = USFM*DELFL
                       FLHAB  = MIN(FLHAB,TEMP)
                       IMATRAA(IP,IS,ID) = SIGN(FLHAB,GTEMP2)/DT4A*JAC
                       LIMFAC            = MIN(ONE,ABS(SIGN(FLHAB,GTEMP2/DT4A))/MAX(THR,ABS(IMATRAA(IP,IS,ID))))
                       IMATDAA(IP,IS,ID) = ZERO ! -FL(IP,ID,IS)
                     ENDIF
                   ENDDO ! ID
                 ENDDO ! IS
                 ACLOC = AC2(IP,:,:)
                 IF (.NOT. LINID) THEN
                   CALL SET_WIND( IP, WIND10, WINDTH )
                   CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
                   CALL SIN_LIN_CAV(IP,WINDTH,FPM,IMATRA,SSINL)
                   IMATRAA(IP,:,:) = IMATRAA(IP,:,:) + SSINL
                 ENDIF
                 IF (ISHALLOW(IP) .EQ. 1) THEN
                   CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)
                   IF (MESTR .GT. 0) THEN
                     CALL TRIADSWAN_NEW (IP, HS, SME01, ACLOC, IMATRA, IMATDA, SSNL3, DSSNL3)
                     DO IS = 1, MSC
                       DO ID = 1, MDC
                         NEWDAC = SSNL3(IS,ID)*DT4A/MAX((1.-DT4A*DSSNL3(IS,ID)),1.)
                         MAXDAC = 0.0081*LIMFAK/(TWO*SPSIG(IS)*WK(IP,IS)**3*CG(IP,IS))*100
                         LIMFAC = ONE/MAX(ONE,NEWDAC/MAXDAC)
                         SC = SIGN(MIN(ABS(NEWDAC),MAXDAC),NEWDAC)/DT4A
                         !SSNL3(IS,ID)  = SC
                         !DSSNL3(IS,ID) = DSSNL3(IS,ID)*LIMFAC
                         !IF (ABS(SC) .GT. THR) WRITE(*,'(2I10,5F20.8)') IS, ID, NEWDAC, MAXDAC, DSSNL3(IS,ID), LIMFAC
                       END DO
                     END DO
                   ENDIF ! MESTR
                   IF (MESBR .GT. 0) CALL SDS_SWB(IP, SME01, KMWAM, ETOT, HS, ACLOC, IMATRA, IMATDA, SSBR, DSSBR)
                   IF (MESBF .GT. 0) CALL SDS_BOTF(IP,ACLOC,IMATRA,IMATDA,SSBF,DSSBF)
                   IMATDAA(IP,:,:) = IMATDAA(IP,:,:) + DSSBR  + DSSNL3 + DSSBF
                   IMATRAA(IP,:,:) = IMATRAA(IP,:,:) + SSBR + SSNL3
                 ENDIF ! ISHALLOW(IP) .EQ. 1
               ENDIF
             ENDIF
           ENDIF
           !IF (IP==TESTNODE) WRITE(*,*) IP, SUM(IMATRAA(IP,:,:)), DEP(IP), IOBP(IP), SUM(SSBR), SUM(SSBF), SUM(SSINL), SUM(SSNL3)
         ENDDO
         !WRITE(*,'(A20,3F15.10)') 'FROM SPECINT', SUM(IMATRAA), SUM(IMATDAA), SUM(AC2)
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
               THWOLD(IP,1) = THWNEW(IP)
               THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))
               U10NEW(IP) = MAX(TWO,SQRT(WINDXY(IP,1)**2+WINDXY(IP,2)**2)) * WINDFAC
               Z0NEW(IP) = Z0OLD(IP,1)
               DO IS = 1, MSC
                 DO ID = 1, MDC
                   FL3(IP,ID,IS) =  AC2(IP,IS,ID) * PI2 * SPSIG(IS)
                   FL(IP,ID,IS)  =  IMATDAA(IP,IS,ID)
                   SL(IP,ID,IS)  =  IMATRAA(IP,IS,ID) * PI2 * SPSIG(IS)
                 END DO
               END DO
               CALL POSTINTRHS (FL3(IP,:,:), FL(IP,:,:), IP, IP, 1, &
     &                          THWOLD(IP,1), USOLD(IP,1), &
     &                          TAUW(IP), Z0OLD(IP,1), &
     &                          ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                          U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                          Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                          SL(IP,:,:), FCONST(IP,:), FMEANWS(IP), MIJ(IP))
               DO IS = 1, MSC
                 DO ID = 1, MDC
                   AC2(IP,IS,ID) = FL3(IP,ID,IS) / PI2 / SPSIG(IS)
                 END DO
               END DO
             END IF !
           ELSE
             IF (LSOUBOUND) THEN ! Source terms on boundary ...
               IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
                 FL = FL3
                 THWOLD(:,1) = THWNEW
                 U10NEW = MAX(TWO,SQRT(WINDXY(:,1)**2+WINDXY(:,2)**2)) * WINDFAC
                 Z0NEW(IP) = Z0OLD(IP,1)
                 THWNEW(IP) = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))
                 DO IS = 1, MSC
                   DO ID = 1, MDC
                     FL3(IP,ID,IS) =  AC2(IP,IS,ID) * PI2 * SPSIG(IS)
                     FL(IP,ID,IS)  =  IMATDAA(IP,IS,ID)
                     SL(IP,ID,IS)  =  IMATRAA(IP,IS,ID) * PI2 * SPSIG(IS)
                   END DO
                 END DO
                 CALL POSTINTRHS (FL3(IP,:,:), FL(IP,:,:), IP, IP, 1, &
     &                            THWOLD(IP,1), USOLD(IP,1), &
     &                            TAUW(IP), Z0OLD(IP,1), &
     &                            ROAIRO(IP,1), ZIDLOLD(IP,1), &
     &                            U10NEW(IP), THWNEW(IP), USNEW(IP), &
     &                            Z0NEW(IP), ROAIRN(IP), ZIDLNEW(IP), &
     &                            SL(IP,:,:), FCONST(IP,:), FMEANWS(IP), MIJ(IP))
                 DO IS = 1, MSC
                   DO ID = 1, MDC
                     AC2(IP,IS,ID) = FL3(IP,ID,IS) / PI2 / SPSIG(IS)
                   END DO
                 END DO
               ENDIF
             ENDIF
           ENDIF
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

         !IF (IP == 1701) WRITE(*,*) SUM(IMATRA), SUM(IMATDA), SUM(ACLOC)

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
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(ZERO,IMATDA(IS,ID)) )
                 IF (LIMIT) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
                 ACLOC(IS,ID) = MAX( ZERO, ACLOC(IS,ID) + NEWDAC )
               END DO
             END DO
             CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE.) ! 2. CALL
             DO IS = 1, MSC_HF(IP)
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(ZERO,IMATDA(IS,ID)) )
                 IF (LIMIT) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
                 ACLOC(IS,ID) = MAX( ZERO, 3./4. * ACOLD(IS,ID) +  1./4. * ACLOC(IS,ID) + 1./4. * NEWDAC)
               END DO
             END DO
             CALL SOURCETERMS(IP, ACLOC, IMATRA, IMATDA, .FALSE.) ! 3. CALL
             DO IS = 1, MSC_HF(IP)
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(ZERO,IMATDA(IS,ID)) )
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
!2do add mean quantities for 
      SUBROUTINE SOURCETERMS (IP, ACLOC, IMATRA, IMATDA, LRECALC)
      USE DATAPOOL, ONLY : MSC, MDC, MNP, WK, LINID, THR, UFRIC
      USE DATAPOOL, ONLY : CD, TAUTOT, TAUWX, TAUWY, AC1, AC2, DEP
      USE DATAPOOL, ONLY : PI2, CG, G9, ZERO, ALPHA_CH, QBLOCAL
      USE DATAPOOL, ONLY : USTDIR, Z0, SMALL, VERYSMALL, MSC_HF
      USE DATAPOOL, ONLY : DDIR, SPSIG, SPDIR, FRINTF, LMAXETOT, LADVTEST
      USE DATAPOOL, ONLY : USTDIR, Z0, SMALL, VERYSMALL, MSC_HF, DDIR
      USE DATAPOOL, ONLY : SPSIG, SPDIR, FRINTF, ONE, RHOA, RHOAW, TAUHF
      USE DATAPOOL, ONLY : TAUW, FR, MESNL, MESIN, MESDS, MESBF, MESBR
      USE DATAPOOL, ONLY : MESTR, ISHALLOW, DS_INCR, IOBP, IOBPD
      USE DATAPOOL, ONLY : LNANINFCHK, DBG, IFRIC, RTIME, DISSIPATION
      USE DATAPOOL, ONLY : AIRMOMENTUM, ONEHALF, NSPEC, RKIND, ISELECT
#ifdef WWM_MPI
      USE DATAPOOL, ONLY : myrank
#endif
         USE SdsBabanin
#ifdef SNL4_TSA
         USE W3SNLXMD
#endif
#ifdef ST41
         USE W3SRC4MD_OLD
#elif ST42
         USE W3SRC4MD
#endif
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP

         REAL(rkind), INTENT(INOUT) :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)

         LOGICAL, INTENT(IN) :: LRECALC

         INTEGER        :: IS, ID, IS0, IK, ITH, IDISP, JU, NZZ
         REAL(rkind)    :: WIND10, WINDTH
         REAL(rkind)    :: FPM
         REAL(rkind)    :: SME01, SME10, KME01, KMWAM, KMWAM2
         REAL(rkind)    :: SME01WS, SME10WS
         REAL(rkind)    :: HS, ETOT, FPMH,FPM4 
         REAL(rkind)    :: LPOW, MPOW, a1, a2, ETOTWS
         REAL(rkind)    :: XRR, XPP, XFLT, XREL, FACP, XFILT
         REAL(rkind)    :: TEMP2(MSC), SSBRL(MSC,MDC)
         REAL(rkind)    :: AWW3(NSPEC), IMATRAWW3(MSC,MDC), IMATDAWW3(MSC,MDC)
         REAL(rkind)    :: EDENS(MSC), KDS(MSC), ABAB(MSC)
         REAL(rkind)    :: WHITECAP(1:4),AKMEAN,XKMEAN,F1MEAN,TMPAC(MDC,MSC),TEMP(MSC), FCONST(MSC)
         REAL(rkind)    :: ACLOC1(MSC,MDC)


         REAL(rkind)    :: SSNL3(MSC,MDC), SSNL4(MSC,MDC), SSINL(MSC,MDC), SSDS(MSC,MDC), DSSNL4(MSC,MDC)
         REAL(rkind)    :: SSBF(MSC,MDC), SSBR(MSC,MDC), SSINE(MSC,MDC), DSSNL3(MSC,MDC), DSSBR(MSC,MDC)
         REAL(rkind)    :: TMP_IN(MSC), TMP_DS(MSC), WN2(MSC*MDC),SPRDD(MDC),AK2VGM1,AKM1, DAM(MSC*MDC)

         REAL(rkind)    :: EMEAN, FMEAN, FMEAN1, WNMEAN, AMAX, AS
         REAL(rkind)    :: FMEANWS, TAUWAX, TAUWAY, XJ, FLLOWEST, GADIAG
         REAL(rkind)    :: IMATDA1D(NSPEC), IMATRA1D(NSPEC), SUMACLOC, IMATRAT(MSC,MDC), BRLAMBDA(NSPEC)
         REAL(rkind)    :: IMATRA_WAM(MDC,MSC), IMATDA_WAM(MDC,MSC), TAILFACTOR, FLOGSPRDM1, SNL3(MSC,MDC), DSNL3(MSC,MDC)
         REAL    :: IMATRA_TSA(MDC,MSC), IMATDA_TSA(MDC,MSC), TMPAC_TSA(MDC,MSC), CG_TSA(MSC), WK_TSA(MSC), DEP_TSA
         REAL    :: XNL(MSC,MDC), DDIAG(MSC,MDC), ACLOC_WRT(MSC,MDC), DEP_WRT, SPSIG_WRT(MSC), SPDIR_WRT(MDC)
         INTEGER :: IERR_WRT

#ifdef TIMINGS 
         REAL(rkind)        :: T1, T2
         REAL(rkind), SAVE  :: TIME1, TIME2, TIME3, TIME4, TIME5, TIME6, TIME7, TIME8, TIME9
#endif

         LOGICAL        :: LWINDSEA(MSC,MDC)
         REAL(rkind)    :: XLCKS(MDC,MSC)

#ifdef TIMINGS
         call MY_WTIME(TIME1)
#endif 

         WIND10 = ZERO 
         SUMACLOC = SUM(ACLOC)

         IF (LMAXETOT .AND. .NOT. LADVTEST .AND. ISHALLOW(IP) .EQ. 1 .AND. .NOT. LRECALC) THEN
           CALL BREAK_LIMIT(IP,ACLOC,SSBRL) ! Miche to reduce stiffness of source terms ...
         END IF

         IDISP = 999

         IF (.NOT. LRECALC) THEN
           ACLOC1=AC1(IP,:,:)
           CALL MEAN_WAVE_PARAMETER(IP,ACLOC1,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2) ! 1st guess ... 
         END IF

         SSINE       = zero
         SSINL       = zero
         SSNL3       = zero
         SSNL4       = zero
         DSSNL4      = zero
         SSBR        = zero
         SSBF        = zero
         IMATRA_WAM  = zero
         IMATDA_WAM  = zero
         TMPAC       = zero
         IMATRA      = zero 
         IMATDA      = zero 
         IMATRAWW3   = zero 
         IMATDAWW3   = zero 
         QBLOCAL(IP) = zero 

#ifdef ST_DEF
         IF (MESDS == 1 .OR. MESIN .EQ. 1) THEN
           DO IS = 1, MSC
             DO ID = 1, MDC
               AWW3(ID + (IS-1) * MDC) = ACLOC(IS,ID) * CG(IP,IS)
             END DO
           END DO
         END IF
         XPP     = 0.15
         XRR     = 0.10
         XFILT  = 0.05
         XPP     = MAX ( 1.E-6_rkind , XPP )
         XRR     = MAX ( 1.E-6_rkind , XRR )
         XREL   = XRR
         XFILT  = MAX ( ZERO , XFILT )
         XFLT   = XFILT
         FACP   = 2*XPP / PI2 * 0.62E-3 * PI2**4 / G9**2

         DO IK=1, NK
           DAM(1+(IK-1)*NTH) = FACP / ( SIG(IK) * WK(IP,IK)**3 )
           WN2(1+(IK-1)*NTH) = WK(IP,IK)
         END DO
!
         DO IK=1, NK
           IS0    = (IK-1)*NTH
           DO ITH=2, NTH
             DAM(ITH+IS0) = DAM(1+IS0)
             WN2(ITH+IS0) = WN2(1+IS0)
           END DO
         END DO
#endif

#ifdef TIMINGS
         call MY_WTIME(TIME2)
#endif 
         IF ((ISELECT .EQ. 1 .OR. ISELECT .EQ. 10 .OR. ISELECT .EQ. 20) .AND. .NOT. LRECALC) THEN
           IF (IOBP(IP) .EQ. 0) THEN
             IF (MESIN == 1) THEN ! Ardhuin et al. 2010
               CALL SET_WIND( IP, WIND10, WINDTH )
               TAUWX(IP) = ZERO
               TAUWY(IP) = ZERO               
               IF (SUMACLOC .LT. THR .AND. WIND10 .GT. THR) THEN
                 CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
                 IF (.NOT. LINID) THEN
                   CALL SIN_LIN_CAV(IP,WINDTH,FPM,IMATRA,SSINL)
                 ENDIF
               ELSE
                 MSC_HF(IP) = MSC
                 AS      = 0. 
#ifdef ST41
                 CALL W3SPR4_OLD ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
                 CALL W3SIN4_OLD ( AWW3, CG(IP,:), WK(IP,:), WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, IMATRA1D, IMATDA1D, LLWS)
                 CALL W3SPR4_OLD ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
                 CALL W3SIN4_OLD ( AWW3, CG(IP,:), WK(IP,:), WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, IMATRA1D, IMATDA1D, LLWS)
#elif ST42
                 CALL W3SPR4 ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
                 CALL W3SIN4 ( IP, AWW3, CG(IP,:), WN2,  WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, IMATRA1D, IMATDA1D, LLWS, BRLAMBDA)
                 CALL W3SPR4 ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)  
                 CALL W3SIN4 ( IP, AWW3, CG(IP,:), WN2, WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, IMATRA1D, IMATDA1D, LLWS, BRLAMBDA) 
#else
                 WRITE(DBG%FHNDL,*) 'NO ST42 or ST41 chosen but MESIN == 1'
                 CALL WWM_ABORT('stop wwm_sourceterms l.169')
#endif
                 CALL ONED2TWOD(IMATRA1D,IMATRAWW3)
                 SSINE = IMATRAWW3
                 CALL ONED2TWOD(IMATDA1D,IMATDAWW3)
                 DO ID = 1, MDC
                   IMATRA(:,ID) = IMATRAWW3(:,ID) / CG(IP,:) + IMATRA(:,ID)
                   IMATDA(:,ID) = IMATDAWW3(:,ID) !/ CG(IP,:) 
                 END DO
               END IF
               IF (LNANINFCHK) THEN
                 IF (SUM(IMATRA) .NE. SUM(IMATRA) .OR. SUM(IMATDA) .NE. SUM(IMATDA)) THEN
                   WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT NORMAL', IP, '   DUE TO SIN', SUM(IMATRA), SUM(IMATDA)
                   CALL WWM_ABORT('NAN AT wwm_sourceterms.F90 l.204')
                 END IF
               ENDIF
             ELSE IF (MESIN == 2) THEN ! Cycle 4, Bidlot et al. ...
               CALL WWM_ABORT('PLEASE USE LSOURCEWAM FOR ECWAM SOURCE TERM FORMULATION') 
             ELSE IF (MESIN == 3) THEN ! Makin & Stam
               CALL SET_WIND( IP, WIND10, WINDTH )
               IFRIC = 4
               CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
               IF (.NOT. LINID) CALL SIN_LIN_CAV( IP, WINDTH, FPM, IMATRA, SSINL)
               CALL SIN_MAKIN( IP, WIND10, WINDTH, KME01,ETOT,ACLOC,IMATRA,IMATDA,SSINE)
             ELSE IF (MESIN == 4) THEN ! Donealan et al.
               CALL SET_WIND( IP, WIND10, WINDTH )
               IFRIC = 1
               CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
               IF (.NOT. LINID) CALL SIN_LIN_CAV( IP, WINDTH,FPM,IMATRA,SSINL)
               CALL SWIND_DBYB (IP,WIND10,WINDTH,IMATRA,SSINE)
             ELSE IF (MESIN == 5) THEN ! Cycle 3
               CALL SET_WIND( IP, WIND10, WINDTH )
               IFRIC = 1
               CALL SET_FRICTION( IP, ACLOC, WIND10, WINDTH, FPM )
               IF (.NOT. LINID) CALL SIN_LIN_CAV( IP, WINDTH,FPM,IMATRA,SSINL)
               CALL SIN_EXP_KOMEN( IP, WINDTH, ACLOC, IMATRA, IMATDA, SSINE )
             END IF ! MESIN
           END IF ! IOBP
           !IF (IDISP == IP) THEN
             !CALL PLOT_SHADED_CONTOUR_POLAR(SPSIG/PI2,SPDIR*RADDEG,MSC,MDC,IMATRA/MAXVAL(IMATRA),50,MSC,MDC,'BEFORE ANY CALL')
             !PAUSE
           !END IF
         END IF ! ISELECT 

#ifdef TIMINGS
         call MY_WTIME(TIME3)
#endif 
         IF (LNANINFCHK) THEN
           IF (SUM(IMATRA) .NE. SUM(IMATRA) .OR. SUM(IMATDA) .NE. SUM(IMATDA)) THEN
             WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, '   DUE TO SIN', SUM(IMATRA), SUM(IMATDA)
             CALL WWM_ABORT('NAN in wwm_sourceterm.F90 l.311')
           END IF
         ENDIF

         IF ((ISELECT.EQ.2 .OR. ISELECT.EQ.10 .OR. ISELECT.EQ.20) .AND. .NOT. LRECALC) THEN
           IF (IOBP(IP) .EQ. 0) THEN
             IF (MESNL .EQ. 1) THEN
               CALL SNL41 (IP, KMWAM, ACLOC, IMATRA, IMATDA, SSNL4, DSSNL4)
             ELSE IF (MESNL .EQ. 2) THEN
               CALL SNL4(IP, KMWAM, ACLOC, IMATRA, IMATDA)
             ELSE IF (MESNL .EQ. 3) THEN
               CALL SNL42(IP, KMWAM, ACLOC, IMATRA, IMATDA)
             ELSE IF (MESNL .EQ. 4) THEN
               CALL SNL43(IP, KMWAM, ACLOC, IMATRA, IMATDA)
             ELSE IF (MESNL .EQ. 5) THEN
               ACLOC_WRT = REAL(ACLOC)
               SPSIG_WRT = REAL(SPSIG)
               SPDIR_WRT = REAL(SPDIR)
               DEP_WRT   = DEP(IP)
               CALL WWMQUAD_WRT (ACLOC_WRT,SPSIG_WRT,SPDIR_WRT,MDC,MSC,DEP_WRT,1,XNL,DDIAG,IERR_WRT)
               IF (IERR_WRT .GT. 0) THEN
                 WRITE (DBG%FHNDL,*) 'XNL_WRT ERROR', IERR_WRT
                 CALL WWM_ABORT('XNL_WRT ERROR')
               ELSE
                 IMATRA(:,:) = IMATRA(:,:) + XNL (:,:)
                 IMATDA(:,:) = IMATDA(:,:) + DDIAG(:,:)
               END IF
             ELSE IF (MESNL .EQ. 6) THEN
               DO IS = 1, MSC
                 DO ID = 1, MDC
                   TMPAC_TSA(ID,IS) = ACLOC(IS,ID) * CG(IP,IS)
                 END DO
               END DO
               CG_TSA = CG(IP,:)
               WK_TSA = WK(IP,:)
               DEP_TSA = DEP(IP)
               NZZ = (MSC*(MSC+1))/2
#ifdef SNL4_TSA
               CALL W3SNLX ( TMPAC_TSA, CG_TSA, WK_TSA, DEP_TSA, NZZ, IMATRA_TSA, IMATDA_TSA)
#endif
               DO IS = 1, MSC
                 DO ID = 1, MDC
                   IMATRA(IS,ID) = IMATRA(IS,ID) + IMATRA_TSA(ID,IS) / CG(IP,IS)
                   IMATDA(IS,ID) = IMATDA(IS,ID) + IMATDA_TSA(ID,IS)
                 END DO
               END DO
             END IF
           END IF
           IF (IDISP == IP) THEN
             !CALL PLOT_SHADED_CONTOUR_POLAR(SPSIG/PI2,SPDIR*RADDEG,MSC,MDC,IMATRA/MAXVAL(IMATRA),50,MSC,MDC,'BEFORE ANY CALL')
             !PAUSE
           END IF
         END IF ! ISELECT

         IF (LNANINFCHK) THEN
           IF (SUM(IMATRA) .NE. SUM(IMATRA) .OR. SUM(IMATDA) .NE. SUM(IMATDA)) THEN
             WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, '   DUE TO SNL4'
             CALL WWM_ABORT('NAN at wwm_sourceterms.F90 l.368')
           END IF
         ENDIF

#ifdef TIMINGS
         call MY_WTIME(TIME4)
#endif 
         IF ((ISELECT.EQ.3 .OR. ISELECT.EQ.10 .OR. ISELECT.EQ.20) .AND. .NOT. LRECALC) THEN

           IMATRAT = IMATRA

           IF (IOBP(IP) .EQ. 0 .OR. IOBP(IP) .EQ. 4) THEN

             IF (MESDS == 1) THEN
#ifdef ST41
               CALL W3SDS4_OLD(AWW3,WK(IP,:),CG(IP,:),UFRIC(IP),USTDIR(IP),DEP(IP),IMATRA1D, IMATDA1D) 
#elif ST42
               CALL W3SDS4(AWW3,WK(IP,:),CG(IP,:),UFRIC(IP),USTDIR(IP),DEP(IP),IMATRA1D,IMATDA1D,BRLAMBDA,WHITECAP)
#endif
               CALL ONED2TWOD(IMATRA1D,IMATRAWW3)
               CALL ONED2TWOD(IMATDA1D,IMATDAWW3)
               DO ID = 1, MDC
                 SSDS(:,ID)   = IMATRAWW3(:,ID) / CG(IP,:)
                 IMATRA(:,ID) = IMATRA(:,ID)+IMATRAWW3(:,ID) / CG(IP,:)
                 IMATDA(:,ID) = IMATDA(:,ID)+IMATDAWW3(:,ID) !/ CG(IP,:)
               END DO
             ELSE IF (MESDS == 2) THEN
               CALL WWM_ABORT('PLEASE USE LSOURCEWAM FOR ECWAM SOURCE TERM FORMULATION')
             ELSE IF (MESDS == 3) THEN
               CALL SDS_NEDWAM_CYCLE4( IP, KMWAM, SME01, ETOT, ACLOC, IMATRA, IMATDA, SSDS  ) ! NEDWAM 
             ELSE IF (MESDS == 4) THEN
               ABAB = 1.
               LPOW = 4.
               MPOW = 4.
               a1  = 0.00000045
               a2  = 0.0000095
!              LPOW = 2.
!              MPOW = 2.
!              a1  = 0.0002
!              a2  = 0.003
!  0.0002 0.003 2.0 2.0 KOM
!  0.00000045 0.0000095 4.0 4.0 BD
               DO IS = 1, MSC
                 EDENS(IS) = 0.
                 DO ID = 1, MDC
                   EDENS(IS) = EDENS(IS) + ACLOC(IS,ID) *  SPSIG(IS) * PI2 * DDIR
                 END DO
               END DO
               CALL CALC_SDS(IP,MSC,EDENS,FR,Kds,ABAB,LPOW,MPOW,a1,a2,ACLOC,IMATRA,IMATDA,SSDS)
!              CALL SSWELL(IP,ETOT,ACLOC,IMATRA,IMATDA,URMSTOP,CG0)
             ELSE IF (MESDS == 5) THEN
               CALL SDS_CYCLE3 ( IP, KMWAM, SME10, ETOT, ACLOC, IMATRA, IMATDA, SSDS ) 
             END IF
           END IF
#ifdef VDISLIN
           IF (IDISP == IP) THEN
             CALL PLOT_SHADED_CONTOUR_POLAR(SPSIG/PI2,SPDIR*RADDEG,MSC,MDC,IMATRAT-IMATRA,50,MSC,MDC,'BEFORE ANY CALL')
           END IF
#endif
         END IF

#ifdef TIMINGS
         call MY_WTIME(TIME5)
#endif 
         IF (((ISELECT.EQ.4 .OR. ISELECT.EQ.10).AND.ISHALLOW(IP).EQ.1) .AND. .NOT. LRECALC) THEN
           IF (SUMACLOC .GT. VERYSMALL) THEN
             IF (MESTR .EQ. 1 ) THEN
               CALL TRIADSWAN_NEW (IP,HS,SME01,ACLOC,IMATRA,IMATDA,SSNL3,DSSNL3)
             ELSE IF (MESTR .EQ. 2) THEN
               CALL SNL31 (IP,HS,SME01,ACLOC,IMATRA,IMATDA,SSNL3)
             ELSE IF (MESTR .EQ. 3) THEN
               CALL SNL32 (IP,HS,SME01,ACLOC,IMATRA,IMATDA,SSNL3)
             ELSE IF (MESTR .EQ. 4) THEN
               CALL SNL33 (IP,HS,SME01,ACLOC,IMATRA,IMATDA,SSNL3)
             ELSE IF (MESTR .EQ. 5) THEN
               CALL TRIADSWAN (IP,HS,SME01,ACLOC,IMATRA,IMATDA,SSNL3)
             ELSE IF (MESTR .EQ. 6) THEN
               CALL WWM_ABORT('NOT READ YET')
               CALL TRIAD_DINGEMANS (IP,ACLOC,IMATRA,IMATDA,SSNL3)
             ELSE IF (MESTR .EQ. 7) THEN
               CALL WWM_ABORT('NOT READ YET')
               CALL snl3ta(ip,snl3,dsnl3)
               SSNL3 = SNL3
               IMATRA = IMATRA + SNL3
               IMATDA = IMATDA + DSNL3
             END IF
           END IF
         END IF

#ifdef TIMINGS
         call MY_WTIME(TIME6)
#endif 
         IF (MESBR .EQ. 1) THEN
           IF (((ISELECT.EQ.5 .OR. ISELECT.EQ.10 .OR. ISELECT.EQ.30) .AND. ISHALLOW(IP) .EQ. 1) .AND. .NOT. LRECALC) THEN
             IF (IOBP(IP) == 0 .OR. IOBP(IP) == 4 .OR. IOBP(IP) == 3) THEN
               CALL SDS_SWB(IP,SME01,KMWAM,ETOT,HS,ACLOC, IMATRA,IMATDA,SSBR,DSSBR)
             END IF
           ENDIF
         END IF
#ifdef TIMINGS
         call MY_WTIME(TIME7)
#endif 
         IF (MESBF .GE. 1) THEN
           IF (((ISELECT.EQ.6 .OR. ISELECT.EQ.10 .OR. ISELECT.EQ.30) .AND. ISHALLOW(IP) .EQ. 1) .AND. .NOT. LRECALC) THEN
             IF (IOBP(IP) == 0 .OR. IOBP(IP) == 4 .OR. IOBP(IP) == 3) THEN
              CALL SDS_BOTF(IP,ACLOC,IMATRA,IMATDA,SSBR,DSSBR)
             END IF
           END IF
         ENDIF
         IF (LNANINFCHK) THEN
           IF (SUM(IMATRA) .NE. SUM(IMATRA) .OR. SUM(IMATDA) .NE. SUM(IMATDA)) THEN
             WRITE(DBG%FHNDL,*) 'NAN AT GRIDPOINT', IP, '   DUE TO SBF' 
             CALL WWM_ABORT('NAN in wwm_sourceterms.F90 at l.481')
           END IF
         ENDIF
#ifdef TIMINGS
        call MY_WTIME(TIME8)
#endif 
!------------------------------------------------------------------------------------------------------------------------!
!-------------------------------- RECALCULATE ALL SOURCE TERMS BASED ON THE NEW SPECTRA ---------------------------------! 
!------------------------------------------------------------------------------------------------------------------------!
         IF (LRECALC .and. IOBP(IP) .EQ. 0) THEN

           DISSIPATION(IP) = 0.
           AIRMOMENTUM(IP) = 0.
           DO ID = 1, MDC
             TMP_DS = ( SSBR(:,ID) + SSBF(:,ID) + SSDS(:,ID) ) * SPSIG * DDIR
             TMP_IN = ( SSINE(:,ID) + SSINL(:,ID) ) * SPSIG * DDIR
             DO IS = 2, MSC
               DISSIPATION(IP) = DISSIPATION(IP) + ONEHALF * ( TMP_DS(IS) + TMP_DS(IS-1) ) * DS_INCR(IS)
               AIRMOMENTUM(IP) = AIRMOMENTUM(IP) + ONEHALF * ( TMP_IN(IS) + TMP_IN(IS-1) ) * DS_INCR(IS)
             END DO
           END DO

           IF (MESIN == 1) THEN
             AS      = 0.
#ifdef ST41
             CALL W3SPR4_OLD ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
             CALL W3SIN4_OLD ( AWW3, CG(IP,:), WK(IP,:), WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, IMATRA1D, IMATDA1D, LLWS)
             CALL W3SPR4_OLD ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
#elif ST42
             CALL W3SPR4 ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
             CALL W3SIN4 ( IP, AWW3, CG(IP,:), WN2,  WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, IMATRA1D, IMATDA1D, LLWS, BRLAMBDA)
             CALL W3SPR4 ( AWW3, CG(IP,:), WK(IP,:), EMEAN, FMEAN, FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS)
#else
             WRITE(DBG%FHNDL,*) 'NO ST42 or ST41 chosen but MESIN == 1'
             CALL WWM_ABORT('stop wwm_sourceterms l.186')
#endif
           ELSEIF (MESIN == 2) THEN
           ELSEIF (MESIN == 3) THEN
           ELSEIF (MESIN == 4) THEN
           ELSEIF (MESIN == 5) THEN
           ENDIF
         ENDIF

#ifdef TIMINGS 
         call MY_WTIME(TIME9)
#ifdef WWM_MPI 
         IF (IP == MNP .AND. myrank == 0 ) THEN
#else if 
         IF (IP == MNP) THEN
#endif
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----SOURCE TIMINGS-----'
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'PREPARATION        ', TIME2-TIME1
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SIN                ', TIME3-TIME2
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SDS                ', TIME4-TIME3
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SNL4               ', TIME5-TIME4
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SNL3               ', TIME6-TIME5
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SBR                ', TIME7-TIME6
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SBF                ', TIME8-TIME7
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'RECALC             ', TIME9-TIME8
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'TOTAL              ', TIME9-TIME1
!           WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '------END-TIMINGS-  ---'
         ENDIF
#endif 

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
