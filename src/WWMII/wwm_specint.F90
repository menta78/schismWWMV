#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SOURCE_INT_EXP()

         USE DATAPOOL
#ifdef ST41
         USE W3SRC4MD_OLD, ONLY : LFIRSTSOURCE
#endif
#ifdef ST42
         USE W3SRC4MD_NEW, ONLY : LFIRSTSOURCE
#endif
         IMPLICIT NONE

         INTEGER        :: IP, IS, ID, NSTEP, MSC_HF
         INTEGER        :: NIT_SIN, NIT_SDS, NIT_SNL4, NIT_SNL3, NIT_SBR, NIT_SBF, NIT_ALL
         REAL(rkind)    :: ACLOC(MSC,MDC), IMATRA(MSC,MDC), IMATDA(MSC,MDC), SSBRL1(MSC,MDC), SSBRL2(MSC,MDC)
         REAL(rkind)    :: WIND10, WINDTH

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(IP,IS,ID,ACLOC) 
         DO IP = 1, MNP 
           IF ((ABS(IOBP(IP)) .NE. 1 .AND. IOBP(IP) .NE. 3)) THEN
             IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
               ACLOC  = AC2(IP,:,:)
               IF (SMETHOD == 1) THEN
!                 DT4S = 0.5 * DT4S 
!                 CALL INT_IP_DYN(IP, 2, .TRUE., DTMIN_SNL4, NDYNITER_SNL4, ACLOC, NIT_SNL4)
!                 DT4S = 2 * DT4S
!                 CALL INT_IP_DYN(IP, 20, .FALSE., DTMIN_DYN , NDYNITER_SIN , ACLOC, NIT_SIN) ! Sbf
!                 DT4S = 0.5 * DT4S
!                 CALL INT_IP_DYN(IP, 2, .TRUE., DTMIN_SNL4, NDYNITER_SNL4, ACLOC, NIT_SNL4)
!                 CALL INT_IP_DYN(IP, 6, .FALSE., DTMIN_SBF , NDYNITER_SBF , ACLOC, NIT_SBF) ! Sbf
!                 CALL INT_IP_DYN(IP, 5, .FALSE., DTMIN_SBR , NDYNITER_SBR , ACLOC, NIT_SBR) ! Sbr
!                 CALL INT_IP_DYN(IP, 4, .FALSE., DTMIN_SNL3, NDYNITER_SNL3, ACLOC, NIT_SNL3)! Snl3 

                 DT4S = 0.5 * DT4S
                 !CALL INT_IP_STAT(IP,2,.TRUE.,ACLOC)
                 CALL RKS_SP3(IP,2,DT4S,.TRUE.,ACLOC)
                 DT4S = 2 * DT4S
                 !CALL INT_IP_STAT(IP,20,.FALSE.,ACLOC)
                 CALL RKS_SP3(IP,20,DT4S,.FALSE.,ACLOC)
                 DT4S = 0.5 * DT4S
                 CALL RKS_SP3(IP,2,DT4S,.TRUE.,ACLOC)
                 !CALL INT_IP_STAT(IP,2,.TRUE.,ACLOC)

                 write(*,*) NIT_SNL4, NIT_SIN, NIT_SBF, NIT_SBR, NIT_SNL3
               ELSE IF (SMETHOD == 2) THEN
                 CALL INT_IP_STAT(IP,10,LLIMT,ACLOC)
               ELSE IF (SMETHOD == 3) THEN
                 CALL RKS_SP3(IP,10,DT4S,LLIMT,ACLOC)
               ELSE IF (SMETHOD == 4) THEN
                 CALL INT_IP_DYN(IP, 10, LLIMT, DTMIN_DYN, NDYNITER, ACLOC, NIT_ALL)
                 write(*,*) NDYNITER, NIT_ALL
               ELSE IF (SMETHOD == 5) THEN ! Full splitting of all source embedded within a dynamic RK-3 Integration ... 
                 CALL INT_IP_DYN(IP, 1, LLIMT, DTMIN_SIN ,   NDYNITER_SIN  , ACLOC, NIT_SIN) ! Sin
                 CALL INT_IP_DYN(IP, 2, LLIMT, DTMIN_SNL4,   NDYNITER_SNL4 , ACLOC, NIT_SNL4)! Snl4b
                 CALL INT_IP_DYN(IP, 3, LLIMT, DTMIN_SDS ,   NDYNITER_SDS  , ACLOC, NIT_SDS) ! Sds
                 CALL INT_IP_DYN(IP, 4, LLIMT, DTMIN_SNL3,   NDYNITER_SNL3 , ACLOC, NIT_SNL3)! Snl3
                 CALL INT_IP_DYN(IP, 5, LLIMT, DTMIN_SBR ,   NDYNITER_SBR  , ACLOC, NIT_SBR) ! Sbr
                 CALL INT_IP_DYN(IP, 6, LLIMT, DTMIN_SBF ,   NDYNITER_SBF  , ACLOC, NIT_SBF) ! Sbf
               END IF
               CALL SOURCETERMS(IP, 1, ACLOC, IMATRA, IMATDA, .TRUE., MSC_HF) ! Update everything based on the new spectrum ...
               IF (LMAXETOT .AND. .NOT. LADVTEST .AND. ISHALLOW(IP) .EQ. 1) THEN
                 CALL BREAK_LIMIT(IP,ACLOC,SSBRL2)
               ENDIF
               AC2(IP,:,:) = ACLOC
             ENDIF
           ELSE !Boundary node ... 
             IF (LSOUBOUND) THEN ! Source terms on boundary ...
               IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
!                 write(*,*) 'calling routine shit bb', ip, iobp(ip)
                 ACLOC  = AC2(IP,:,:)
                 IF (SMETHOD == 1) THEN
                   CALL INT_IP_DYN(IP, 6, .FALSE., DTMIN_SBF , NDYNITER_SBF , ACLOC, NIT_SBF) ! Sbf
                   CALL INT_IP_DYN(IP, 5, .FALSE., DTMIN_SBR , NDYNITER_SBR , ACLOC, NIT_SBR) ! Sbr
                   CALL INT_IP_DYN(IP, 4, .FALSE., DTMIN_SNL3, NDYNITER_SNL3, ACLOC, NIT_SNL3)! Sbf
                   write(*,*) NIT_SBF, NIT_SBR, NIT_SNL3
                 ELSE IF (SMETHOD == 2) THEN
                   CALL INT_IP_STAT(IP,10,LLIMT,ACLOC)
                 ELSE IF (SMETHOD == 3) THEN
                   CALL RKS_SP3(IP,10,DT4S,LLIMT,ACLOC)
                 ELSE IF (SMETHOD == 4) THEN
                   CALL INT_IP_DYN(IP, 10, LLIMT, DTMIN_DYN,NDYNITER, ACLOC, NIT_ALL)
                 ELSE IF (SMETHOD == 5) THEN ! Full splitting of all source embedded within a dynamic RK-3 Integration ... 
                   CALL INT_IP_DYN(IP, 1, LLIMT, DTMIN_SIN ,   NDYNITER_SIN  , ACLOC, NIT_SIN) ! Sin
                   CALL INT_IP_DYN(IP, 2, LLIMT, DTMIN_SNL4,   NDYNITER_SNL4 , ACLOC, NIT_SNL4)! Snl4b
                   CALL INT_IP_DYN(IP, 3, LLIMT, DTMIN_SDS ,   NDYNITER_SDS  , ACLOC, NIT_SDS) ! Sds
                   CALL INT_IP_DYN(IP, 4, LLIMT, DTMIN_SNL3,   NDYNITER_SNL3 , ACLOC, NIT_SNL3)! Snl3
                   CALL INT_IP_DYN(IP, 5, LLIMT, DTMIN_SBR ,   NDYNITER_SBR  , ACLOC, NIT_SBR) ! Sbr
                   CALL INT_IP_DYN(IP, 6, LLIMT, DTMIN_SBF ,   NDYNITER_SBF  , ACLOC, NIT_SBF) ! Sbf
                 END IF
                 CALL SOURCETERMS(IP, 1, ACLOC, IMATRA, IMATDA, .TRUE., MSC_HF) ! Update everything based on the new spectrum ...
                 IF (LMAXETOT .AND. .NOT. LADVTEST .AND. ISHALLOW(IP) .EQ. 1) THEN
                   CALL BREAK_LIMIT(IP,ACLOC,SSBRL2)
                   AC2(IP,:,:) = ACLOC
                 END IF
                 AC2(IP,:,:) = ACLOC
               ENDIF
             ELSE
               IF (LMAXETOT .AND. .NOT. LADVTEST .AND. ISHALLOW(IP) .EQ. 1) THEN ! limit wave height on the boundary ...
                 CALL BREAK_LIMIT(IP,ACLOC,SSBRL2)
                 AC2(IP,:,:) = ACLOC
               ENDIF
             ENDIF 
           ENDIF
         ENDDO
#if defined ST41 || defined ST42
         LFIRSTSOURCE = .FALSE.
#endif
         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SOURCE_INT_IMP()

         USE DATAPOOL
#ifdef ST41
         USE W3SRC4MD_OLD, ONLY : LFIRSTSOURCE
#endif
#ifdef ST42
         USE W3SRC4MD_NEW, ONLY : LFIRSTSOURCE
#endif
         IMPLICIT NONE

         INTEGER :: IP, ID, ISELECT, MSC_HF

         REAL(rkind)    :: ACLOC(MSC,MDC)
         REAL(rkind)    :: IMATRA(MSC,MDC), IMATDA(MSC,MDC), SSBRL(MSC,MDC)

!$OMP WORKSHARE
         IMATDAA = 0.
         IMATRAA = 0.
!$OMP END WORKSHARE

         ISELECT = 10

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(IP,ACLOC,IMATDA,IMATRA)
         DO IP = 1, MNP
           IF ((ABS(IOBP(IP)) .NE. 1 .OR. IOBP(IP) .NE. 3) .AND. LSOUBOUND) THEN
             IF ( DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 2) THEN
               MSC_HF = MSC
               ACLOC = AC2(IP,:,:)
               CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA,.FALSE.,MSC_HF) 
               IMATDAA(IP,:,:) = IMATDA
               IMATRAA(IP,:,:) = IMATRA
             END IF !
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
      SUBROUTINE INT_IP_STAT(IP,ISELECT,LIMITER,ACLOC)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP, ISELECT
         LOGICAL, INTENT(IN) :: LIMITER
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)
         REAL(rkind)                :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER :: IS, ID, MSC_HF

         REAL(rkind)    :: NPF(MSC)
         REAL(rkind)    :: OLDAC
         REAL(rkind)    :: NEWDAC, NEWAC(MSC,MDC)
         REAL(rkind)    :: MAXDAC, CONST, SND

         CONST = PI2**2*3.0*1.0E-7*DT4S*SPSIG(MSC)
         SND   = PI2*5.6*1.0E-3

         IF (LIMITER) THEN
           IF (MELIM == 1) THEN
             NPF = 0.0081_rkind*LIMFAK/(2._rkind*SPSIG*WK(IP,:)**3.0_rkind*CG(IP,:))
           ELSE IF (MELIM == 2) THEN
             NPF = LIMFAK*ABS((CONST*(MAX(UFRIC(IP),G9*SND/SPSIG)))/(SPSIG**3*WK(IP,:)))
           END IF
         END IF

         CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA, .FALSE., MSC_HF)  ! 1. CALL

!         if (SUM(ACLOC) .NE. SUM(ACLOC)) STOP 'NAN l. 174 wwm_specint.F90'

         DO IS = 1, MSC
           MAXDAC = NPF(IS)
           DO ID = 1, MDC
             OLDAC  = ACLOC(IS,ID)
             NEWDAC = IMATRA(IS,ID) * DT4S / ( 1.0 - DT4S * MIN(ZERO,IMATDA(IS,ID))) 
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
      SUBROUTINE INT_IP_ECMWF(IP,ISELECT,LIMITER,ACLOC)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP, ISELECT
         LOGICAL, INTENT(IN) :: LIMITER
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)
         REAL(rkind)                :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)
         REAL(rkind)                :: IMATRA_WAM(MSC,MDC), IMATDA_WAM(MSC,MDC)


!         DELT = DT4S
!         XIMP = 1.0
!         DELT5 = XIMP*DELT

!         DO IS=1,MSC
!           DELFL(IS) = COFRM4(IS)*DELT
!         ENDDO

!         USFM = USNEW*MAX(FMEANWS,FMEAN)

!         DO IS=1,MSC
!           TEMP(IS) = USFM*DELFL(IS)
!         ENDDO

!         IF(ISHALLO.EQ.1) THEN
!           DO IS=1,MSC
!             TEMP2(IS) = FRM5(IS)
!           ENDDO
!         ELSE
!           DO IS=1,MSC
!             AKM1      = ONE/WK(IP,IS)
!             AK2VGM1   = AKM1**2/CG(IP,IS)
!             TEMP2(IS) = AKM1*AK2VGM1
!           ENDDO
!         ENDIF

!         DO ID=1,MDC
!           DO IS=1,MSC 
!             IMATRA_WAM(IS,ID) = IMATRA(IS,ID) * PI2 * SPSIG(IS)
!             IMATDA_WAM(IS,ID) = IMATDA(IS,ID) * PI2 * SPSIG(IS)
!             GTEMP1 = MAX((1.-DELT5*IMATDA_WAM(IS,ID)),1.)
!             GTEMP2 = DELT*IMATRA_WAM(IS,ID)/GTEMP1
!             FLHAB = ABS(GTEMP2)
!             FLHAB = MIN(FLHAB,TEMP(IS))
!             ACLOC(IP,IS,ID) = ACLOC(IS,ID) + SIGN(FLHAB,GTEMP2)
!             FLLOWEST = VERYSMALL 
!             ACLOC(IS,ID) = MAX(ACLOC(IS,ID),FLLOWEST)
!           ENDDO
!         ENDDO

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE RKS_SP3(IP,ISELECT,DTSII,LIMITER,ACLOC)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP, ISELECT
         LOGICAL, INTENT(IN) :: LIMITER
         REAL(rkind), INTENT(IN)    :: DTSII
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)

         INTEGER :: IS, ID, MSC_HF
         REAL(rkind)    :: ACOLD(MSC,MDC)
         REAL(rkind)    :: NPF(MSC)
         REAL(rkind)    :: NEWDAC, MAXDAC, CONST, SND
         REAL(rkind)    :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         CONST = PI2**2*3.0*1.0E-7*DTSII*SPSIG(MSC)
         SND   = PI2*5.6*1.0E-3
 
         IMATRA = 0.
         IMATDA = 0.

         IF (LIMITER) THEN
           IF (MELIM == 1) THEN
             NPF = 0.0081_rkind*LIMFAK/(TWO*SPSIG*WK(IP,:)**THREE*CG(IP,:))
           ELSE IF (MELIM == 2) THEN
             NPF = LIMFAK*ABS((CONST*(MAX(UFRIC(IP),G9*SND/SPSIG)))/(SPSIG**3*WK(IP,:)))
           END IF
         END IF

         CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA, .FALSE., MSC_HF)  ! 1. CALL

         ACOLD = ACLOC

         DO IS = 1, MSC
           MAXDAC = NPF(IS)
           DO ID = 1, MDC
             NEWDAC = IMATRA(IS,ID) * DTSII / ( ONE - DTSII * MIN(0._rkind,IMATDA(IS,ID)) )
             IF (LIMITER) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
             ACLOC(IS,ID) = MAX( ZERO, ACLOC(IS,ID) + NEWDAC )
           END DO
         END DO

         !WRITE(*,*) '1 RK-TVD', SUM(ACOLD), SUM(ACLOC)

         CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA, .FALSE., MSC_HF) ! 2. CALL

         DO IS = 1, MSC
           MAXDAC = NPF(IS)
           DO ID = 1, MDC
             NEWDAC = IMATRA(IS,ID) * DTSII / ( ONE - DTSII * MIN(ZERO,IMATDA(IS,ID)) )
             IF (LIMITER) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC) 
             ACLOC(IS,ID) = MAX( ZERO, 0.75_rkind * ACOLD(IS,ID) +  0.25_rkind * ACLOC(IS,ID) + 0.25_rkind * NEWDAC)
           END DO
         END DO

         !WRITE(*,*) '2 RK-TVD', SUM(ACOLD), SUM(ACLOC)

         CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA, .FALSE., MSC_HF) ! 3. CALL

         DO IS = 1, MSC
           MAXDAC = NPF(IS)
           DO ID = 1, MDC
             NEWDAC = IMATRA(IS,ID) * DTSII / ( 1.0 - DTSII * MIN(ZERO,IMATDA(IS,ID)) )
             IF (LIMITER) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
             ACLOC(IS,ID) = MAX( ZERO, ONE/THREE * ACOLD(IS,ID) + TWO/THREE * ACLOC(IS,ID) + TWO/THREE * NEWDAC)
           END DO
         END DO

         !WRITE(*,*) '3 RK-TVD', SUM(ACOLD), SUM(ACLOC)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INT_IP_DYN(IP, ISELECT, LIMIT, DTMIN, ITRMX, ACLOC, ITER)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)        :: IP, ISELECT, ITRMX
         INTEGER, INTENT(OUT)       :: ITER
         LOGICAL,INTENT(IN)          :: LIMIT
         REAL(rkind), INTENT(IN)    :: DTMIN
         REAL(rkind), INTENT(INOUT) :: ACLOC(MSC,MDC)
         REAL(rkind)                :: IMATRA(MSC,MDC), IMATDA(MSC,MDC)

         INTEGER :: IS, ID, MSC_HF

         REAL(rkind)    :: ACOLD(MSC,MDC)
         REAL(rkind)    :: NPF(MSC)
         REAL(rkind)    :: TMP1, TMP2, TMP3
         REAL(rkind)    :: NEWDAC, CONST, SND
         REAL(rkind)    :: MAXDAC, DTMAX, DTTOT, DTLEFT, DT4SI

         REAL(rkind), PARAMETER :: MAXDTFAC = VERYLARGE 

         LOGICAL :: LSTABLE

         CONST = PI2**2*3.0*1.0E-7*DT4S*SPSIG(MSC)
         SND   = PI2*5.6*1.0E-3

         IF (MELIM == 1) THEN
           NPF = 0.0081_rkind*LIMFAK/(two*SPSIG*WK(IP,:)**3*CG(IP,:))
         ELSE IF (MELIM == 2) THEN
           NPF = LIMFAK*ABS((CONST*(MAX(UFRIC(IP),G9*SND/SPSIG)))/(SPSIG**3*WK(IP,:)))
         END IF

         DTTOT = 0.
         ITER  = 0

         DO WHILE ( DTTOT < DT4S )

           ACOLD = ACLOC
           IF (ITER == 0) DT4SI = DT4S
           IF (ITER == 0) DTLEFT = DT4S
           ITER  = ITER + 1
           DTMAX =  DT4S 

           CALL SOURCETERMS(IP, ISELECT,ACLOC, IMATRA, IMATDA, .FALSE., MSC_HF)  ! 1. CALL
        
           DO ID = 1, MDC
             DO IS = 1, MSC_HF
               IF (ABS(IMATRA(IS,ID)) .GT. VERYSMALL) THEN                
                 DTMAX = MIN(DTMAX,MIN(DT4S,NPF(IS)/ABS(IMATRA(IS,ID))))
               END IF
             END DO
           END DO

           DT4SI  = DTMAX
           DT4SI  = MAX(DTMIN,DT4SI) ! This is not entirely stable !!!
           DTLEFT = DT4S - DTTOT

           IF ( DTLEFT > THR .AND. DTLEFT < DT4SI) THEN
             DT4SI = (DT4S - DTTOT)
           ELSE IF ( DTLEFT .GE. DT4SI .AND. ITER .EQ. ITRMX) THEN
             DT4SI = DTLEFT 
             LSTABLE = .FALSE. 
           END IF 

           DTTOT = DTTOT + DT4SI

           IF ( DT4SI .LT. DTMAX ) THEN
             LSTABLE = .TRUE. 
           ELSE
             LSTABLE = .FALSE.
           END IF

           IF (LSTABLE) THEN
             DO IS = 1, MSC_HF
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(ZERO,IMATDA(IS,ID)) )
                 ACLOC(IS,ID) = MAX( ZERO, ACLOC(IS,ID) + NEWDAC )
               END DO
             END DO
             CALL SOURCETERMS(IP,ISELECT,ACLOC,IMATRA,IMATDA,.FALSE.,MSC_HF) ! 2. CALL
             DO IS = 1, MSC_HF
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(ZERO,IMATDA(IS,ID)) )
                 ACLOC(IS,ID) = MAX( ZERO, 3._rkind/4._rkind * ACOLD(IS,ID) + 1._rkind/4._rkind * ACLOC(IS,ID) + 1._rkind/4._rkind * NEWDAC)
               END DO
             END DO
             CALL SOURCETERMS(IP,ISELECT,ACLOC,IMATRA,IMATDA, .FALSE., MSC_HF) ! 3. CALL
             DO IS = 1, MSC_HF
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(ZERO,IMATDA(IS,ID)) )
                 ACLOC(IS,ID) = MAX( ZERO,  1._rkind/3._rkind * ACOLD(IS,ID) + 2._rkind/3._rkind * ACLOC(IS,ID) + 2._rkind/3._rkind * NEWDAC)
               END DO
             END DO
           ELSE ! .NOT. LSTABLE
             DO IS = 1, MSC_HF
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(0._rkind,IMATDA(IS,ID)) )
                 IF (LIMIT) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
                 ACLOC(IS,ID) = MAX( ZERO, ACLOC(IS,ID) + NEWDAC )
               END DO
             END DO
             CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA, .FALSE., MSC_HF) ! 2. CALL
             DO IS = 1, MSC_HF
               MAXDAC = NPF(IS)
               DO ID = 1, MDC
                 NEWDAC = IMATRA(IS,ID) * DT4SI / ( 1.0 - DT4SI * MIN(0._rkind,IMATDA(IS,ID)) )
                 IF (LIMIT) NEWDAC = SIGN(MIN(MAXDAC,ABS(NEWDAC)),NEWDAC)
                 ACLOC(IS,ID) = MAX( ZERO, 3./4. * ACOLD(IS,ID) +  1./4. * ACLOC(IS,ID) + 1./4. * NEWDAC)
               END DO
             END DO
             CALL SOURCETERMS(IP, ISELECT, ACLOC, IMATRA, IMATDA, .FALSE., MSC_HF) ! 3. CALL
             DO IS = 1, MSC_HF
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
