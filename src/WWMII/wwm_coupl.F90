#include "wwm_functions.h"
#undef DEBUG_WWM
#define DEBUG_WWM
#if !defined PGMCL_COUPLING
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_PIPES_TIMOR()
         USE DATAPOOL
         IMPLICIT NONE
!
! open pipe data files for coupling
!
           LSEWL = .TRUE.
           LSECU = .TRUE.
           WRITE(*,'("+TRACE...",A)') 'OPEN PIPE'
!          Pipes that are read by the wave model
           OPEN(1000,file='p_velx.dat'  ,form='unformatted', action='read')
           OPEN(1001,file='p_vely.dat'  ,form='unformatted', action='read')
           OPEN(1002,file='p_lev.dat'   ,form='unformatted', action='read')
           OPEN(1003,file='p_bot.dat'   ,form='unformatted', action='read')
           OPEN(1004,file='p_time.dat'  ,form='unformatted', action='write')
           OPEN(1005,file='p_dthyd.dat',form='unformatted', action='read')
!          Pipes that are written by the wave modell
           OPEN(101 ,file='p_stressx.dat' ,form='unformatted', action='write')
           OPEN(102 ,file='p_stressy.dat' ,form='unformatted', action='write')
           OPEN(103 ,file='p_waveh.dat'   ,form='unformatted', action='write')
           OPEN(104 ,file='p_wavet.dat'   ,form='unformatted', action='write')
           OPEN(105 ,file='p_waved.dat'   ,form='unformatted', action='write')
           OPEN(106 ,file='p_wavekm.dat'  ,form='unformatted', action='write')
           OPEN(107 ,file='p_wavetp.dat'  ,form='unformatted', action='write')
           OPEN(108 ,file='p_wavekp.dat'  ,form='unformatted', action='write')
           OPEN(109 ,file='p_orbit.dat'   ,form='unformatted', action='write')
           OPEN(110 ,file='p_stokesx.dat' ,form='unformatted', action='write')
           OPEN(111 ,file='p_stokesy.dat' ,form='unformatted', action='write')
           OPEN(112 ,file='p_windx.dat'   ,form='unformatted', action='write')
           OPEN(113 ,file='p_windy.dat'   ,form='unformatted', action='write')
!          Pipes writen as ergzus.bin for XFN
           OPEN(2003, file='pipe_wave.bin' ,form='unformatted' ,status = 'unknown')
           OPEN(2004, file='pipe_orbi.bin' ,form='unformatted' ,status = 'unknown')
           OPEN(2005, file='pipe_stok.bin' ,form='unformatted' ,status = 'unknown')
           WRITE(*,'("+TRACE...",A)') 'END OPEN PIPE'
!
! write coupling tyme step
!
           WRITE(1004) MAIN%DTCOUP
           CALL FLUSH(1004)
!
! read coupling tyme step
!
           READ(1005)  MAIN%DTCUR

           IF(ABS(MAIN%DTCOUP-INT(MAIN%DTCOUP/MAIN%DTCUR)*MAIN%DTCUR).GT. .001) THEN
             write(*,*) 'TIME STEP OF THE hydraulic flow MODEL CANNOT BE DIVIDIED WITHOUT A REST'
             write(*,*)'dt Stroemung (s) =',MAIN%DTCUR, ',  dt Kopplung (s) = ',MAIN%DTCOUP
!             CALL WWM_ABORT('dt Stroemungsmodell muss gerades Vielfaches des Kopplungs-dt sein')
           END IF

           WRITE(*,'("+TRACE... DTCUR and DTCOUP",A)') MAIN%DTCUR, MAIN%DTCOUP

       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE TERMINATE_PIPES_TIMOR()
      USE DATAPOOL
      IMPLICIT NONE
      close(1000)
      close(1001)
      close(1002)
      close(1003)
      close(1004)
      close(1005)
      close(101)
      close(102)
      close(103)
      close(104)
      close(105)
      close(106)
      close(107)
      close(108)
      close(109)
      close(110)
      close(111)
      close(112)
      close(113)
      close(2003)
      close(2004)
      close(2005)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PIPE_TIMOR_IN(K)
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)  :: K

         IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN
           WRITE(*,'("+TRACE...",A)') 'READING PIPE'
           WRITE(*,'("+TRACE...",A)') 'END READING PIPE'
         END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PIPE_TIMOR_OUT(K)
         IMPLICIT NONE

         INTEGER, INTENT(IN)  :: K

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_PIPES_SHYFEM()

         USE DATAPOOL
         IMPLICIT NONE
!
! open pipe data files for coupling
!
           LSEWL = .TRUE.
           LSECU = .TRUE.

           WRITE(*,'("+TRACE...",A)') 'OPEN PIPE'
!          Pipes that are read by the wave model
           OPEN(1000,file='p_velx.dat'    ,form='unformatted', action='read')
           OPEN(1001,file='p_vely.dat'    ,form='unformatted', action='read')
           OPEN(1002,file='p_lev.dat'     ,form='unformatted', action='read')
           OPEN(1003,file='p_bot.dat'     ,form='unformatted', action='read')
           OPEN(1004,file='p_zeta3d.dat'  ,form='unformatted', action='read')
!          Pipes that are written by the wave model
           OPEN(101 ,file='p_stressx.dat' ,form='unformatted', action='write')
           OPEN(102 ,file='p_stressy.dat' ,form='unformatted', action='write')
           OPEN(142 ,file='p_stresxy.dat' ,form='unformatted', action='write')  !ccf
           OPEN(103 ,file='p_waveh.dat'   ,form='unformatted', action='write')
           OPEN(104 ,file='p_wavet.dat'   ,form='unformatted', action='write')
           OPEN(105 ,file='p_waved.dat'   ,form='unformatted', action='write')
           OPEN(106 ,file='p_wavekm.dat'  ,form='unformatted', action='write')
           OPEN(107 ,file='p_wavetp.dat'  ,form='unformatted', action='write')
           OPEN(108 ,file='p_wavekp.dat'  ,form='unformatted', action='write')
           OPEN(109 ,file='p_orbit.dat'   ,form='unformatted', action='write')
           OPEN(110 ,file='p_stokesx.dat' ,form='unformatted', action='write')
           OPEN(111 ,file='p_stokesy.dat' ,form='unformatted', action='write')

           IF (.NOT. LWINDFROMWWM) THEN
!            *** WIND FROM SHYFEM *** ccf
             OPEN(112,file='p_windx.dat'  ,form='unformatted', action='read')
             OPEN(113,file='p_windy.dat'  ,form='unformatted', action='read')
            ELSE
!            *** WIND FROM WWM *** ccf
             OPEN(112 ,file='p_windx.dat',form='unformatted', action='write')
             OPEN(113 ,file='p_windy.dat',form='unformatted', action='write')
           END IF

           OPEN(114 ,file='p_cd.dat' ,form='unformatted', action='write')
           OPEN(115 ,file='p_jpress.dat' ,form='unformatted', action='write')

           WRITE(*,'("+TRACE...",A)') 'END OPEN PIPE'

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE TERMINATE_PIPES_SHYFEM()
      USE DATAPOOL
      IMPLICIT NONE
      close(1000)
      close(1001)
      close(1002)
      close(1003)
      close(1004)
      close(101)
      close(102)
      close(142)
      close(103)
      close(104)
      close(105)
      close(106)
      close(107)
      close(108)
      close(109)
      close(110)
      close(111)
      close(112)
      close(113)
      close(114)
      close(115)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PIPE_SHYFEM_IN(K)
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)  :: K
         INTEGER              :: IP, IL

         IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN

           WRITE(STAT%FHNDL,'("+TRACE...",A)') 'READING PIPE'
           IF (.NOT. LWINDFROMWWM) THEN
!      *** WIND FROM SHYFEM *** ccf
             DO IP = 1, MNP
               READ(1000) CURTXY(IP,1)
               READ(1001) CURTXY(IP,2)
               READ(1002) WATLEV(IP)
               READ(1003) WLDEP(IP)
               READ(1003) NLEV(IP)
               DO IL = 1,NLVT
                 READ(1004) SHYFZETA(IL,IP)
               END DO
               READ(112) WINDXY(IP,1)
               READ(113) WINDXY(IP,2)
             END DO
            ELSE
!      *** WIND FROM WWM *** ccf
             DO IP = 1, MNP
               READ(1000) CURTXY(IP,1)
               READ(1001) CURTXY(IP,2)
               READ(1002) WATLEV(IP)
               READ(1003) WLDEP(IP)
               READ(1003) NLEV(IP)
               DO IL = 1,NLVT
                 READ(1004) SHYFZETA(IL,IP)
               END DO
             END DO
             WRITE(STAT%FHNDL,*) 'CHECK MAX UX,UY,H'
             WRITE(STAT%FHNDL,*) MAXVAL(CURTXY(:,1)), MAXVAL(CURTXY(:,2)), MAXVAL(WATLEV)
             WRITE(STAT%FHNDL,*) 'CHECK MIN UX,UY,H'
             WRITE(STAT%FHNDL,*) MINVAL(CURTXY(:,1)), MINVAL(CURTXY(:,2)), MINVAL(WATLEV)
             WRITE(2001,*) 'CHECK MAX UX,UY,H'
             WRITE(2001,*) MAXVAL(CURTXY(:,1)), MAXVAL(CURTXY(:,2)), MAXVAL(WATLEV)
             WRITE(2001,*) 'CHECK MIN UX,UY,H'
             WRITE(2001,*) MINVAL(CURTXY(:,1)), MINVAL(CURTXY(:,2)), MINVAL(WATLEV)
           END IF
           WRITE(STAT%FHNDL,'("+TRACE...",A)') 'END READ PIPE WWM'
           LCALC = .TRUE.

         END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PIPE_SHYFEM_OUT(K)
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)  :: K

         INTEGER         :: IL, IP
         REAL(rkind)     :: ACLOC(MSC,MDC)
         REAL(rkind)     :: HS,WLM,LPP
         REAL(rkind)     :: SME01,SME10,KME01,KMWAM,KMWAM2,URSELL,UBOT
         REAL(rkind)     :: TMBOT,ETOT,FP,CP,KPP,DM,DSPR,PEAKDSPR,PEAKDM
         REAL(rkind)     :: ORBITAL, USTOKES, USTOKES_X, USTOKES_Y
         REAL(rkind)     :: BOTEXPER, ETOTS, ETOTC, DPEAK
         REAL(rkind)     :: FPP, TPP, CPP, WNPP, CGPP, TPPD,KPPD,CGPD,CPPD 

           IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN

             WRITE(STAT%FHNDL,'("+TRACE...",A)') 'WRITING PIPE'
!
!            *** COMPUTE RADIATION STRESSES 2D OR 3D *** ccf
!
             IF (LCPL) THEN
               CALL RADIATION_STRESS_SHYFEM
             END IF

             IF (.NOT. LWINDFROMWWM) THEN
               DO IP = 1, MNP
                 DO IL = 1, NLVT
                   WRITE(101)  SXX3D(IL,IP)             !ccf
                   WRITE(102)  SYY3D(IL,IP)             !ccf
                   WRITE(142)  SXY3D(IL,IP)             !ccf
                 END DO
                 CALL FLUSH(101)
                 CALL FLUSH(102)
                 CALL FLUSH(142)                        !ccf
                 ACLOC = AC2(IP,:,:)
!2do check ustokes ...
                 CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)
                 CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT)
                 CALL MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,MSC,ETOTS,ETOTC,DM,DSPR)
                 CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)

                 WRITE(103) HS
                 CALL FLUSH(103)
                 WRITE(104)  SME01
                 CALL FLUSH(104)
                 WRITE(105)  DM
                 CALL FLUSH(105)
                 WRITE(106)  KMWAM
                 CALL FLUSH(106)
                 WRITE(107)  TPP
                 CALL FLUSH(107)
                 WRITE(108)  KPP
                 CALL FLUSH(108)
                 WRITE(109)  ORBITAL
                 CALL FLUSH(109)
!AR: Delete this stuff 
                 WRITE(110)  USTOKES_X
                 CALL FLUSH(110)
                 WRITE(111)  USTOKES_Y
                 CALL FLUSH(111)
               END DO
             ELSE
               DO IP = 1, MNP
                 DO IL = 1, NLVT
                   WRITE(101)  SXX3D(IL,IP)             !ccf
                   WRITE(102)  SYY3D(IL,IP)             !ccf
                   WRITE(142)  SXY3D(IL,IP)             !ccf
                 END DO
                 CALL FLUSH(101)
                 CALL FLUSH(102)
                 CALL FLUSH(142)                        !ccf
                 CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)
                 CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT)
                 CALL MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,MSC,ETOTS,ETOTC,DM,DSPR)
                 CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
                 WRITE(103) HS
                 CALL FLUSH(103)
                 WRITE(104)  SME01
                 CALL FLUSH(104)
                 WRITE(105)  DM
                 CALL FLUSH(105)
                 WRITE(106)  KMWAM
                 CALL FLUSH(106)
                 WRITE(107)  TPP
                 CALL FLUSH(107)
                 WRITE(108)  KPP
                 CALL FLUSH(108)
                 WRITE(109)  ORBITAL
                 CALL FLUSH(109)
!AR: Delete this stuff 
                 !WRITE(110)  USTOKES_X
                 CALL FLUSH(110)
                 !WRITE(111)  USTOKES_Y
                 CALL FLUSH(111)
                 WRITE(112)  WINDXY(IP,1)
                 CALL FLUSH(112)
                 WRITE(113)  WINDXY(IP,2)
                 CALL FLUSH(113)
                 WRITE(114)  CD(IP) 
                 CALL FLUSH(114)
               END DO
             END IF
             WRITE(STAT%FHNDL,'("+TRACE...",A)') 'END WRITING PIPE'
          END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_PIPES_ROMS()
      USE DATAPOOL
      IMPLICIT NONE
!
! open pipe data files for coupling
!
      LSEWL = .TRUE.
      LSECU = .TRUE.
      WRITE(*,'("+TRACE...",A)') 'OPEN PIPE ROMS'
!     Pipes that are read by the wave model
      OPEN(1000,file='pipe/ExchRW'  ,form='unformatted', action='read')
      WRITE(DBG%FHNDL,*) 'WWM: open pipe ExchImport'
!     Pipes that are written by the wave modell
      OPEN(101 ,file='pipe/ExchWR' ,form='unformatted', action='write')
      WRITE(DBG%FHNDL,*) 'WWM: open pipe ExchExport'
      WRITE(*,'("+TRACE...",A)') 'END OPEN PIPE ROMS'

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE TERMINATE_PIPES_ROMS()
      USE DATAPOOL
      IMPLICIT NONE
      close(1000)
      close(101)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PIPE_ROMS_IN(K,IFILE,IT)
         USE DATAPOOL
# ifdef WWM_MPI
         USE ELFE_GLBL, ONLY : iplg, np_global
         USE elfe_msgp, only : myrank, nproc, comm, istatus, ierr
# endif
         IMPLICIT NONE
# ifdef WWM_MPI
         include 'mpif.h'
# endif
         INTEGER, INTENT(IN)  :: K,IFILE,IT
         INTEGER              :: IP
# ifdef WWM_MPI
         REAL(rkind), allocatable :: WINDXY_TOT(:,:), CURTXY_TOT(:,:), WATLEV_TOT(:)
         real(rkind), allocatable :: rbuf_real(:)
         integer idx, iProc
# endif
         IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN
           WRITE(*,'("+TRACE...",A)') 'READING PIPE'
# ifndef WWM_MPI
           DO IP = 1, MNP
             READ(1000) WINDXY(IP,1), WINDXY(IP,2), CURTXY(IP,1), CURTXY(IP,2), WATLEV(IP)
           END DO
# else
           allocate(WINDXY_TOT(np_global,2))
           allocate(CURTXY_TOT(np_global,2))
           allocate(WATLEV_TOT(np_global))
           allocate(rbuf_real(np_global*5))
           IF (myrank.eq.0) THEN
             DO IP = 1, np_global
               READ(1000) WINDXY_TOT(IP,1), WINDXY_TOT(IP,2), CURTXY_TOT(IP,1), CURTXY_TOT(IP,2), WATLEV_TOT(IP)
             END DO
             DO IP=1,np_global
               idx=idx+1
               rbuf_real(idx)=WINDXY_TOT(IP,1)
               idx=idx+1
               rbuf_real(idx)=WINDXY_TOT(IP,2)
               idx=idx+1
               rbuf_real(idx)=CURTXY_TOT(IP,1)
               idx=idx+1
               rbuf_real(idx)=CURTXY_TOT(IP,2)
               idx=idx+1
               rbuf_real(idx)=WATLEV_TOT(IP)
             END DO
             DO iProc=2,nproc
               CALL MPI_SEND(rbuf_real,np_global*5,MPI_REAL8, iProc-1, 196, comm, ierr)
             END DO
           ELSE
             CALL MPI_RECV(rbuf_real,np_global*5,MPI_REAL8, 0, 196, comm, istatus, ierr)
             idx=0
             DO IP=1,np_global
               idx=idx+1
               WINDXY_TOT(IP,1)=rbuf_real(idx)
               idx=idx+1
               WINDXY_TOT(IP,2)=rbuf_real(idx)
               idx=idx+1
               CURTXY_TOT(IP,1)=rbuf_real(idx)
               idx=idx+1
               CURTXY_TOT(IP,2)=rbuf_real(idx)
               idx=idx+1
               WATLEV_TOT(IP)=rbuf_real(idx)
             END DO
           END IF
           DO IP = 1, MNP
             WINDXY(IP,:)=WINDXY_TOT(iplg(IP),:)
             CURTXY(IP,:)=CURTXY_TOT(iplg(IP),:)
             WATLEV(IP)=WATLEV_TOT(iplg(IP))
           END DO
           deallocate(rbuf_real)
           deallocate(WINDXY_TOT)
           deallocate(CURTXY_TOT)
           deallocate(WATLEV_TOT)
# endif
           WRITE(*,'("+TRACE...",A)') 'END READING PIPE'
         END IF
!
!2do make a initialization section for ROMS and WWM
!
         IF (K == 1) CALL INITIAL_CONDITION(IFILE,IT)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PIPE_ROMS_OUT(K)
         USE DATAPOOL
# ifdef WWM_MPI
         USE ELFE_GLBL, ONLY : iplg, np_global
         USE elfe_msgp, only : myrank, comm, ierr, rtype
# endif
         IMPLICIT NONE
# ifdef WWM_MPI
         include 'mpif.h'
# endif
!
         INTEGER, INTENT(IN)  :: K 

         INTEGER              :: IP
         REAL(rkind)          :: ACLOC(MSC,MDC)
         REAL(rkind)          :: HS,WLM,LPP,FPP,CPP,BOTEXPER
         REAL(rkind)          :: URSELL,UBOT,TM,TM01,TM10
         REAL(rkind)          :: TMBOT,FP,CP,KPP,DM,DSPR,ORBITAL,ETOTS,ETOTC,WNPP,TPP,CGPP
         REAL(rkind)          :: PEAKDSPR, PEAKDM, KME, SME, HSWE, HSLIM, TM02, KLM, DPEAK
         REAL(rkind)          :: TPPD,KPPD,CGPD,CPPD
# ifdef WWM_MPI
         REAL(rkind), allocatable :: OUTT(:,:), OUTT_TOT(:,:)
         REAL(rkind)    :: TP, MS, DEG
# endif
         IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN
# ifndef WWM_MPI
           DO IP = 1, MNP
             ACLOC = AC2(IP,:,:)
             CALL MEAN_PARAMETER(IP,ACLOC,MSC,HS,TM01,TM02,TM10,KLM,WLM)
             CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT)
             CALL MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,MSC,ETOTS,ETOTC,DM,DSPR)
             CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
! HS, HSWE, HSLIM  ! - Significant wave height (m) -- HS
! DM  ! - Wave direction (degrees)
! TPP ! - Surface wave relative peak period (s) -- TP
! WLM, KME, SME ! - Average Wave Length [m] - LME
! ORBITAL(IP) ! - Wave bottom orbital velocity (m/s)
! TMBOT ! - Bottom wave period (s)
! DISSIPATION(IP) ! - Wave energy dissipation (W/m2)
! QBLOCAL(IP) ! - Percent of breakig waves (nondimensional)
! DSPR ! - directional spreading
! PEAKDSPR ! - peak directional spreading
! PEAKDM ! - Peak direction
!
!AR: what for you need this HSWE and HSLIM and so on ... what is exactly SME in your definition ...
!AR: I have deleted them ...
!
             HSWE=0
             HSLIM=0
             WRITE(101)  HS, HSWE,                                     &
     &                   HSLIM, DM,                                    &
     &                   TPP, WLM,                                     &
     &                   KLM, TM01,                                    &
     &                   ORBITAL, TMBOT,                               &
     &                   DISSIPATION(IP), QBLOCAL(IP),                 &
     &                   DSPR, PEAKDSPR,                               &
     &                   PEAKDM, TM02
!             WRITE(101)  HS, DM, TPP, WLM, KLM, TM01, ORBITAL, TMBOT, &
!     &                   DISSIPATION(IP), QBLOCAL(IP), DSPR, PEAKDSPR,&
!     &                   PEAKDM, TM02
             CALL FLUSH(101)
           END DO
# else
           allocate(OUTT(np_global,16))
           allocate(OUTT_TOT(np_global,16))
           OUTT=0
           DO IP = 1, MNP
             ACLOC = AC2(IP,:,:)
             CALL MEAN_PARAMETER(IP,ACLOC,MSC,HS,TM01,TM02,TM10,KLM,WLM)
             CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT)
             CALL MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,MSC,ETOTS,ETOTC,DM,DSPR)
             CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
             HSWE=0
             HSLIM=0
             OUTT(iplg(IP), 1)=HS
             OUTT(iplg(IP), 2)=HSWE
             OUTT(iplg(IP), 3)=HSLIM
             OUTT(iplg(IP), 4)=DM
             OUTT(iplg(IP), 5)=TPP
             OUTT(iplg(IP), 6)=WLM
             OUTT(iplg(IP), 7)=KLM
             OUTT(iplg(IP), 8)=TM01
             OUTT(iplg(IP), 9)=ORBITAL
             OUTT(iplg(IP),10)=TMBOT
             OUTT(iplg(IP),11)=DISSIPATION(IP)
             OUTT(iplg(IP),12)=QBLOCAL(IP)
             OUTT(iplg(IP),13)=DSPR
             OUTT(iplg(IP),14)=PEAKDSPR
             OUTT(iplg(IP),15)=PEAKDM
             OUTT(iplg(IP),16)=TM02
           END DO
           call mpi_reduce(OUTT,OUTT_TOT,NP_GLOBAL*16,rtype,MPI_SUM,0,comm,ierr)
           IF (myrank.eq.0) THEN
             DO IP=1,NP_GLOBAL
               OUTT_TOT(IP,:)=OUTT_TOT(IP,:)/nwild_gb(IP)
               WRITE(101) OUTT_TOT(IP, 1), OUTT_TOT(IP, 2),              &
     &                    OUTT_TOT(IP, 3), OUTT_TOT(IP, 4),              &
     &                    OUTT_TOT(IP, 5), OUTT_TOT(IP, 6),              &
     &                    OUTT_TOT(IP, 7), OUTT_TOT(IP, 8),              &
     &                    OUTT_TOT(IP, 9), OUTT_TOT(IP,10),              &
     &                    OUTT_TOT(IP,11), OUTT_TOT(IP,12),              &
     &                    OUTT_TOT(IP,13), OUTT_TOT(IP,14),              &
     &                    OUTT_TOT(IP,15), OUTT_TOT(IP,16)
             END DO
           END IF
           deallocate(OUTT)
           deallocate(OUTT_TOT)
# endif
         END IF
         WRITE(*,*) 'export WWM: ending of writing data'
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef PGMCL_COUPLING
# ifdef WWM_MPI
      SUBROUTINE WWM_CreateMatrixPartition
        USE DATAPOOL
        USE elfe_glbl, only : np_global, iplg
        USE elfe_msgp, only : myrank, istatus,itype, ierr
        USE mod_coupler
        implicit none
        integer, allocatable :: rbuf_int(:)
        integer, allocatable :: TheIndex(:)
        integer, allocatable :: NumberNode(:), NumberTrig(:)
        integer, allocatable :: All_LocalToGlobal(:,:)
        integer i, eIdx, iProc, MNPloc, MNEloc, idx
        integer IPc, IP
#  ifdef DEBUG_WWM
        integer MinValIndex, MinValIndexInv, eVal
#  endif
        allocate(MatrixBelongingWAV(np_global, NnodesWAV))
        allocate(NumberNode(NnodesWAV))
        allocate(NumberTrig(NnodesWAV))
        IF (myrank.ne.MyRankLocal) THEN
          CALL WWM_ABORT('die from ignominious death')
        END IF
        IF (MyRankLocal.eq.0) THEN
          allocate(All_LocalToGlobal(np_global, NnodesWAV))
          All_LocalToGlobal=0
          MatrixBelongingWAV=0
          DO i=1,MNP
            eIdx=iplg(i)
            MatrixBelongingWAV(eIdx,1)=1
            All_LocalToGlobal(i,1)=eIdx
          ENDDO
          NumberNode(1)=MNP
          DO iProc=2,NnodesWAV
            allocate(rbuf_int(2))
            CALL MPI_RECV(rbuf_int,2,itype, iProc-1, 194, WAV_COMM_WORLD, istatus, ierr)
            MNPloc=rbuf_int(1)
            MNEloc=rbuf_int(2)
            NumberNode(iProc)=MNPloc
            NumberTrig(iProc)=MNEloc
            deallocate(rbuf_int)
!
            allocate(rbuf_int(MNPloc))
            CALL MPI_RECV(rbuf_int,MNPloc,itype, iProc-1, 195, WAV_COMM_WORLD, istatus, ierr)
            DO IP=1,MNPloc
              eIdx=rbuf_int(IP)
              MatrixBelongingWAV(eIdx,iProc)=1
              All_LocalToGlobal(IP,iProc)=eIdx
            END DO
            deallocate(rbuf_int)
          END DO
!
          allocate(rbuf_int(np_global*NnodesWAV))
          idx=0
          DO iProc=1,NnodesWAV
            DO IP=1,np_global
              idx=idx+1
              rbuf_int(idx)=MatrixBelongingWAV(IP,iProc)
            END DO
          END DO
          DO iProc=2,NnodesWAV
            CALL MPI_SEND(rbuf_int,np_global*NnodesWAV,itype, iProc-1, 196, WAV_COMM_WORLD, ierr)
          END DO
          deallocate(rbuf_int)
        ELSE
          allocate(rbuf_int(2))
          rbuf_int(1)=MNP
          rbuf_int(2)=MNE
          CALL MPI_SEND(rbuf_int,2,itype, 0, 194, WAV_COMM_WORLD, ierr)
          deallocate(rbuf_int)

          allocate(rbuf_int(MNP))
          DO i=1,MNP
            rbuf_int(i)=iplg(i)
          END DO
          CALL MPI_SEND(rbuf_int,MNP,itype, 0, 195, WAV_COMM_WORLD, ierr)
          deallocate(rbuf_int)
!
          allocate(rbuf_int(np_global*NnodesWAV))
          CALL MPI_RECV(rbuf_int,np_global*NnodesWAV,itype, 0, 196, WAV_COMM_WORLD, istatus, ierr)
          idx=0
          DO iProc=1,NnodesWAV
            DO i=1,np_global
              idx=idx+1
              MatrixBelongingWAV(i,iProc)=rbuf_int(idx)
            END DO
          END DO
          deallocate(rbuf_int)
        ENDIF

        allocate(ReindexPerm_wav(MNP))
        allocate(ReindexPermInv_wav(MNP))
        allocate(TheIndex(np_global))
        TheIndex=0
        DO IP=1,MNP
          IPc=iplg(IP)
          TheIndex(IPc)=IP
        END DO
        idx=0
        DO IPc=1,np_global
          IP=TheIndex(IPc)
          IF (IP.gt.0) THEN
            idx=idx+1
            ReindexPerm_wav(idx)=IP
            ReindexPermInv_wav(IP)=idx
          END IF
        END DO
#  ifdef DEBUG_WWM
        MinValIndex=300
        MinValIndexInv=300
        DO IP=1,MNP
          eVal=ReindexPerm_wav(IP)
          IF (eVal.lt.MinValIndex) THEN
            MinValIndex=eVal
          END IF
          eVal=ReindexPermInv_wav(IP)
          IF (eVal.lt.MinValIndexInv) THEN
            MinValIndexInv=eVal
          END IF
        END DO
        WRITE(DBG%FHNDL,*) 'MinValIndex(Dir,Inv)=', MinValIndex, MinValIndexInv
#  endif
        deallocate(TheIndex)
        deallocate(NumberNode)
        deallocate(NumberTrig)
      END SUBROUTINE
# else
      SUBROUTINE WWM_CreateMatrixPartition
        USE DATAPOOL
        USE mod_coupler
        IMPLICIT NONE
        integer IP
        allocate(MatrixBelongingWAV(MNP, 1))
        allocate(ReindexPerm_wav(MNP))
        allocate(ReindexPermInv_wav(MNP))
        DO IP=1,MNP
          ReindexPerm_wav(IP)=IP
          ReindexPermInv_wav(IP)=IP
        ENDDO
        MatrixBelongingWAV=1
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WWM_DeallocMatrixPartition
      USE DATAPOOL
      USE mod_coupler
      IMPLICIT NONE
      deallocate(MatrixBelongingWAV)
      deallocate(ReindexPerm_wav)
      deallocate(ReindexPermInv_wav)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
# ifdef WWM_MPI
      SUBROUTINE WWM_CreateGlobalLON_LAT_ListTrig
        USE mod_coupler
        USE elfe_msgp, only : parallel_abort
        USE elfe_glbl, only : np_global, ne_global
        implicit none
        integer i, j, k, iegb, stat
        character(len=256) :: RHEADER
        integer nb1, nb2
        integer TheId
        allocate(LONtrig_wav(np_global))
        allocate(LATtrig_wav(np_global))
        allocate(ListTrig_wav(ne_global,3))
        TheId=13000
        open(TheId,file='hgrid.gr3',status='old',iostat=stat)
        read(TheId,*) RHEADER
        read(TheId,*) nb1, nb2
        IF ((nb1.ne.ne_global).or.(nb2.ne.np_global)) THEN
          WRITE(DBG%FHNDL,*) 'nb1=', nb1, ' ne_global=', ne_global
          WRITE(DBG%FHNDL,*) 'nb2=', nb2, ' np_global=', np_global
          CALL WWM_ABORT('Inconsistency')
        END IF
        do i=1,np_global
          read(TheId,*) j, LONtrig_wav(i), LATtrig_wav(i)
        enddo
        do i=1,ne_global
          read(TheId,*) iegb,j,(ListTrig_wav(iegb,k),k=1,3)
          IF (j/=3) then
            call wwm_abort('WWM_CreateGlobalLON_LAT_ListTrig error')
          END IF
        enddo
        close(TheId)
      END SUBROUTINE
# else
      SUBROUTINE WWM_CreateGlobalLON_LAT_ListTrig
        USE mod_coupler
        USE DATAPOOL
        implicit none
        allocate(LONtrig_wav(MNP))
        allocate(LATtrig_wav(MNP))
        allocate(ListTrig_wav(MNE,3))
        LONtrig_wav=XP
        LATtrig_wav=YP
        ListTrig_wav=INE
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ROMS_COUPL_INITIALIZE
        USE DATAPOOL, only : MNP, rkind, DEP, XP, YP, np_total, ne_total
        USE mod_coupler
        USE PGMCL_LIBRARY
        USE pgmcl_interp
        USE elfe_msgp, only : istatus, ierr, itype, myrank
        implicit none
        logical DoNearest
        integer, allocatable :: rbuf_int(:)
        integer IP, iNodeSel, idx, eRankRecv
        real(rkind) eDiff, AbsDiff, SumDep1, SumDep2, SumDiff
        real(rkind) minBathy, maxBathy
!        character(len=40) :: FileSave1, FileSave2
!        character(len=3) :: eStrFi
        real(rkind) SumDepReceive
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 1, rnk=', myrank
# endif
        CALL SetComputationalNodes(ArrLocal, NnodesWAV, OCNid)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 1.2, rnk=', myrank
# endif
        CALL WWM_CreateMatrixPartition
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 1.3, rnk=', myrank
# endif
        CALL WWM_CreateGlobalLON_LAT_ListTrig
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 2, rnk=', myrank
# endif
        IF (MyRankLocal.eq.0) THEN
          CALL M2M_send_fem_r8(ArrLocal, OCNid,                         &
     &        np_total, ne_total, LONtrig_wav, LATtrig_wav,             &
     &        ListTrig_wav)
          CALL M2M_send_node_partition(ArrLocal, OCNid,                 &
     &        np_total, NnodesWAV, MatrixBelongingWAV)
        ENDIF
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 3, rnk=', myrank
# endif
        CALL M2M_recv_r8_grid(ArrLocal, OCNid,                          &
     &   xi_rho, eta_rho, LON_rho_ocn, LAT_rho_ocn, MSK_rho_ocn)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 4, rnk=', myrank
# endif
        CALL M2M_recv_r8_grid(ArrLocal, OCNid,                          &
     &   xi_u, eta_u, LON_u_ocn, LAT_u_ocn, MSK_u_ocn)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 5, rnk=', myrank
# endif
        CALL M2M_recv_r8_grid(ArrLocal, OCNid,                          &
     &   xi_v, eta_v, LON_v_ocn, LAT_v_ocn, MSK_v_ocn)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 6, rnk=', myrank
# endif
        CALL M2M_recv_node_partition(ArrLocal, OCNid,                   &
     &   NnodeRho, NnodesOCN, MatrixBelongingOCN_rho)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 7, rnk=', myrank
# endif
        CALL M2M_recv_node_partition(ArrLocal, OCNid,                   &
     &   NnodeU, NnodesOCN, MatrixBelongingOCN_u)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 8, rnk=', myrank
# endif
        CALL M2M_recv_node_partition(ArrLocal, OCNid,                   &
     &   NnodeV, NnodesOCN, MatrixBelongingOCN_v)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 9, rnk=', myrank
# endif
        allocate(rbuf_int(1))
        eRankRecv=ArrLocal % ListFirstRank(OCNid)
        CALL MPI_RECV(rbuf_int,1,itype, eRankRecv, 103, MPI_COMM_WORLD, istatus, ierr)
        Nlevel=rbuf_int(1)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 10, rnk=', myrank
# endif
        ALLOCATE(z_w_loc(0:Nlevel))
        ALLOCATE(eUSTOKES_loc(Nlevel))
        ALLOCATE(eVSTOKES_loc(Nlevel))
        deallocate(rbuf_int)
        DoNearest=.TRUE.
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 11, rnk=', myrank
# endif
        CALL GetString(MyRankGlobal, eStr)
        FileSave_OCNtoWAV_rho='InterpSave_OCNtoWAV_rho'
        CALL SAVE_CreateInterpolationSparseMatrix_r8_FD_2_r8_FE(        &
     &    FileSave_OCNtoWAV_rho, eStr, mMat_OCNtoWAV_rho, DoNearest,    &
     &    xi_rho, eta_rho, LON_rho_ocn, LAT_rho_ocn, MSK_rho_ocn,       &
     &    ne_total, ListTrig_wav,                                       &
     &    np_total, LONtrig_wav, LATtrig_wav)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 12, rnk=', myrank
# endif
        FileSave_OCNtoWAV_u='InterpSave_OCNtoWAV_u'
        CALL SAVE_CreateInterpolationSparseMatrix_r8_FD_2_r8_FE(        &
     &    FileSave_OCNtoWAV_u, eStr, mMat_OCNtoWAV_u, DoNearest,        &
     &    xi_u, eta_u, LON_u_ocn, LAT_u_ocn, MSK_u_ocn,                 &
     &    ne_total, ListTrig_wav,                                       &
     &    np_total, LONtrig_wav, LATtrig_wav)
        FileSave_OCNtoWAV_v='InterpSave_OCNtoWAV_v'
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 13, rnk=', myrank
# endif
        CALL SAVE_CreateInterpolationSparseMatrix_r8_FD_2_r8_FE(        &
     &    FileSave_OCNtoWAV_v, eStr, mMat_OCNtoWAV_v, DoNearest,        &
     &    xi_v, eta_v, LON_v_ocn, LAT_v_ocn, MSK_v_ocn,                 &
     &    ne_total, ListTrig_wav,                                       &
     &    np_total, LONtrig_wav, LATtrig_wav)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 14  aaa, rnk=', myrank
# endif
        CALL MPI_INTERP_GetSystemOutputSide(ArrLocal, OCNid, WAVid,     &
     &    MatrixBelongingOCN_rho, MatrixBelongingWAV,                   &
     &    mMat_OCNtoWAV_rho, TheArr_OCNtoWAV_rho)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'nbNeedTot=', TheArr_OCNtoWAV_rho % nbNeedTot
        WRITE(DBG%FHNDL,*) 'nbProc=', TheArr_OCNtoWAV_rho % nbProc
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, WAV, step 14, rnk=', myrank
# endif
        CALL MPI_INTERP_GetSystemOutputSide(ArrLocal, OCNid, WAVid,     &
     &    MatrixBelongingOCN_u, MatrixBelongingWAV,                     &
     &    mMat_OCNtoWAV_u, TheArr_OCNtoWAV_u)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 15, rnk=', myrank
# endif
        CALL MPI_INTERP_GetSystemOutputSide(ArrLocal, OCNid, WAVid,     &
     &    MatrixBelongingOCN_v, MatrixBelongingWAV,                     &
     &    mMat_OCNtoWAV_v, TheArr_OCNtoWAV_v)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 16, rnk=', myrank
# endif
        CALL DeallocSparseMatrix(mMat_OCNtoWAV_rho)
        CALL DeallocSparseMatrix(mMat_OCNtoWAV_u)
        CALL DeallocSparseMatrix(mMat_OCNtoWAV_v)
        FileSave_WAVtoOCN_rho='InterpSave_WAVtoOCN_rho'
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 17, rnk=', myrank
# endif
        CALL SAVE_CreateInterpolationSparseMatrix_r8_FE_2_r8_FD(        &
     &    FileSave_WAVtoOCN_rho, eStr, mMat_WAVtoOCN_rho, DoNearest,    &
     &    ne_total, ListTrig_wav,                                       &
     &    np_total, LONtrig_wav, LATtrig_wav,                           &
     &    xi_rho, eta_rho, LON_rho_ocn, LAT_rho_ocn, MSK_rho_ocn)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 18, rnk=', myrank
# endif
        FileSave_WAVtoOCN_u='InterpSave_WAVtoOCN_u'
        CALL SAVE_CreateInterpolationSparseMatrix_r8_FE_2_r8_FD(        &
     &    FileSave_WAVtoOCN_u, eStr, mMat_WAVtoOCN_u, DoNearest,        &
     &    ne_total, ListTrig_wav,                                       &
     &    np_total, LONtrig_wav, LATtrig_wav,                           &
     &    xi_u, eta_u, LON_u_ocn, LAT_u_ocn, MSK_u_ocn)
        FileSave_WAVtoOCN_v='InterpSave_WAVtoOCN_v'
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 19, rnk=', myrank
# endif
        CALL SAVE_CreateInterpolationSparseMatrix_r8_FE_2_r8_FD(        &
     &    FileSave_WAVtoOCN_v, eStr, mMat_WAVtoOCN_v, DoNearest,        &
     &    ne_total, ListTrig_wav,                                       &
     &    np_total, LONtrig_wav, LATtrig_wav,                           &
     &    xi_v, eta_v, LON_v_ocn, LAT_v_ocn, MSK_v_ocn)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 20, rnk=', myrank
# endif
        CALL MPI_INTERP_GetSystemInputSide(ArrLocal, WAVid, OCNid,      &
     &    MatrixBelongingWAV, MatrixBelongingOCN_rho,                   &
     &    mMat_WAVtoOCN_rho, TheArr_WAVtoOCN_rho)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 21, rnk=', myrank
# endif
        CALL MPI_INTERP_GetSystemInputSide(ArrLocal, WAVid, OCNid,      &
     &    MatrixBelongingWAV, MatrixBelongingOCN_u,                     &
     &    mMat_WAVtoOCN_u, TheArr_WAVtoOCN_u)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 22, rnk=', myrank
# endif
        CALL MPI_INTERP_GetSystemInputSide(ArrLocal, WAVid, OCNid,      &
     &    MatrixBelongingWAV, MatrixBelongingOCN_v,                     &
     &    mMat_WAVtoOCN_v, TheArr_WAVtoOCN_v)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 23, rnk=', myrank
# endif
        CALL DeallocSparseMatrix(mMat_WAVtoOCN_rho)
        CALL DeallocSparseMatrix(mMat_WAVtoOCN_u)
        CALL DeallocSparseMatrix(mMat_WAVtoOCN_v)
        allocate(ang_rho(MNP))
        allocate(dep_rho(MNP))
        allocate(A_wav_rho_3D(MNP, Nlevel+1))
        allocate(A_wav_rho(MNP))
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 24, rnk=', myrank
# endif
# ifdef DUMMY_COUPLING
        CALL MPI_INTERP_Dummy_Recv(ArrLocal, 203, OCNid)
# else
        CALL MPI_INTERP_RECV_r8(TheArr_OCNtoWAV_rho, 203, A_wav_rho)
        DO idx=1,MNP
          IP=ReindexPerm_wav(idx)
          ang_rho(IP)=A_wav_rho(idx)
        END DO
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 25, rnk=', myrank
# endif
# ifdef DUMMY_COUPLING
        CALL MPI_INTERP_Dummy_Recv(ArrLocal, 217, OCNid)
# else
        CALL MPI_INTERP_RECV_r8(TheArr_OCNtoWAV_rho, 217, A_wav_rho)
# endif
# ifdef DEBUG_WWM
        SumDepReceive=0
# endif
        DO idx=1,MNP
          IP=ReindexPerm_wav(idx)
          dep_rho(IP)=A_wav_rho(idx)
# ifdef DEBUG_WWM
          SumDepReceive=SumDepReceive + abs(A_wav_rho(idx))
# endif
        END DO
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'SumDepReceive=', SumDepReceive
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, WAV, step 33, rnk=', myrank
        AbsDiff=0
        SumDep1=0
        SumDep2=0
        SumDiff=0
        minBathy=140000
        maxBathy=0
        DO IP=1,MNP
          IF (dep_rho(IP).lt.minBathy) THEN
            minBathy=dep_rho(IP)
          END IF
          IF (dep_rho(IP).gt.maxBathy) THEN
            maxBathy=dep_rho(IP)
          END IF
        END DO
        WRITE(DBG%FHNDL,*) 'dep_rho, min=', minBathy, ' max=', maxBathy
        minBathy=140000
        maxBathy=0
        DO IP=1,MNP
          IF (DEP(IP).lt.minBathy) THEN
            minBathy=DEP(IP)
          END IF
          IF (DEP(IP).gt.maxBathy) THEN
            maxBathy=DEP(IP)
          END IF
        END DO
        WRITE(DBG%FHNDL,*) 'DEP, min=', minBathy, ' max=', maxBathy


        iNodeSel=-1
!        CALL MyGetString(MyRankGlobal, eStrFi)
!        FileSave1='DEP_infos' // eStrFi
!        FileSave2='Lookup_infos' // eStrFi
!        open(745, FILE=TRIM(FileSave1))
!        open(746, FILE=TRIM(FileSave2))
        DO IP=1,MNP
          DEP(IP)=dep_rho(IP)
          eDiff=abs(dep_rho(IP) - DEP(IP))
          SumDiff=SumDiff + eDiff
          SumDep1=SumDep1 + dep_rho(IP)
          SumDep2=SumDep2 + DEP(IP)
          IF (eDiff.gt.AbsDiff) THEN
            AbsDiff=eDiff
            iNodeSel=IP
          END IF
          IF ((DEP(IP).ge.200).and.(eDiff.ge.10)) THEN
            WRITE(DBG%FHNDL,*) 'AD, IP=', IP, dep_rho(IP), DEP(IP)
            WRITE(DBG%FHNDL,*) 'AD, xp, yp=', XP(IP), YP(IP)
          END IF
        END DO
!        close(745)
!        close(746)
!        WRITE(DBG%FHNDL,*) 'AD, AbsDiff=', AbsDiff
!        WRITE(DBG%FHNDL,*) 'AD, IP=', iNodeSel, dep_rho(iNodeSel), DEP(iNodeSel)
!        WRITE(DBG%FHNDL,*) 'AD, xp, yp=', XP(iNodeSel), YP(iNodeSel)
        WRITE(DBG%FHNDL,*) 'AD, SumDep1=', SumDep1, ' SumDep2=', SumDep2
        WRITE(DBG%FHNDL,*) 'AD, SumDiff=', SumDiff
# endif
        allocate(z_r(Nlevel))
        allocate(PartialU1(Nlevel))
        allocate(PartialV1(Nlevel))
        allocate(PartialU2(Nlevel))
        allocate(PartialV2(Nlevel))
        allocate(z_w_wav(MNP, 0:Nlevel))
        allocate(U_wav(MNP, Nlevel))
        allocate(V_wav(MNP, Nlevel))
        allocate(USTOKES_wav(MNP, Nlevel))
        allocate(VSTOKES_wav(MNP, Nlevel))
        allocate(ZETA_CORR(MNP))
        allocate(J_PRESSURE(MNP))
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'End ROMS_COUPL_INITIALIZE'
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ROMS_COUPL_DEALLOCATE
        USE DATAPOOL, only : MNP, rkind, DEP, XP, YP, np_total, ne_total
        USE mod_coupler
        USE PGMCL_LIBRARY
        USE pgmcl_interp
        USE elfe_msgp, only : istatus, ierr, itype, myrank
        implicit none
        logical DoNearest
        integer, allocatable :: rbuf_int(:)
        integer IP, iNodeSel, idx, eRankRecv
        real(rkind) eDiff, AbsDiff, SumDep1, SumDep2, SumDiff
        real(rkind) minBathy, maxBathy
        real(rkind) SumDepReceive
        CALL SetComputationalNodes(ArrLocal, NnodesWAV, OCNid)
        CALL WWM_DeallocMatrixPartition
        deallocate(LONtrig_wav, LATtrig_wav, ListTrig_wav)
        deallocate(LON_rho_ocn, LAT_rho_ocn, MSK_rho_ocn)
        deallocate(LON_u_ocn, LAT_u_ocn, MSK_u_ocn)
        deallocate(LON_v_ocn, LAT_v_ocn, MSK_v_ocn)
        deallocate(MatrixBelongingOCN_rho)
        deallocate(MatrixBelongingOCN_u)
        deallocate(MatrixBelongingOCN_v)
        DEALLOCATE(z_w_loc)
        DEALLOCATE(eUSTOKES_loc)
        DEALLOCATE(eVSTOKES_loc)
        CALL DEALLOCATE_Arr(TheArr_OCNtoWAV_rho)
        CALL DEALLOCATE_Arr(TheArr_OCNtoWAV_u)
        CALL DEALLOCATE_Arr(TheArr_OCNtoWAV_v)
        CALL DEALLOCATE_Arr(TheArr_WAVtoOCN_rho)
        CALL DEALLOCATE_Arr(TheArr_WAVtoOCN_u)
        CALL DEALLOCATE_Arr(TheArr_WAVtoOCN_v)
        deallocate(ang_rho)
        deallocate(dep_rho)
        deallocate(A_wav_rho_3D)
        deallocate(A_wav_rho)
        deallocate(z_r)
        deallocate(PartialU1)
        deallocate(PartialV1)
        deallocate(PartialU2)
        deallocate(PartialV2)
        deallocate(z_w_wav)
        deallocate(U_wav)
        deallocate(V_wav)
        deallocate(USTOKES_wav)
        deallocate(VSTOKES_wav)
        deallocate(ZETA_CORR)
        deallocate(J_PRESSURE)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE STOKES_STRESS_INTEGRAL_ROMS
        USE mod_coupler
        USE DATAPOOL
        implicit none
        integer IP, k, ID, IS
        real(rkind) eF1, eF2, eDelta, TheInt, eDep, eHeight
        real(rkind) eFrac, eFracB, eQuot, TheIntChk
        real(rkind) eFct, eQuot1, eQuot2, eQuot3, eScal, eZeta
        real(rkind) eOmega, eMult, MFACT, kD, eSinc
        real(rkind) USTOKES1, USTOKES2, USTOKES3
        real(rkind) VSTOKES1, VSTOKES2, VSTOKES3
        real(rkind) USTOKESpart, VSTOKESpart, eJPress
        real(rkind) ACLOC, eWk, eSigma, eLoc, eSinhkd, eSinh2kd, eSinhkd2
        real(rkind) zMid, HS, ETOT, MinVal_MFACT, MaxVal_MFACT, MaxVal_eQuot1
        real(rkind) MinVal_MFACT_gl, MaxVal_MFACT_gl, MaxHS, SumHS, AvgHS
        real(rkind) WLM, KLM, AvgStokesNormA, AvgStokesNormB
        real(rkind) PPTAIL, CETAIL, CKTAIL
        real(rkind) ETOT1, EKTOT
        real(rkind) eQuotDispersion, eMaxAC, TotSumAC, eQuotAC, eQuotK
        real(rkind) StokesNormA, StokesNormB, cPhase
        integer IDsel, ISsel, SelectedK
        logical DoTail
        real(rkind) SumNormStokesA(Nlevel), SumNormStokesB(Nlevel)
        real(rkind) SumZetaCorr, MaxZetaCorr, AvgZetaCorr
        real(rkind) eMinMfact, eMaxMfact, SelectedHS
        real(rkind) MaxStokesNorm, MaxValSinc, StokesNorm, SelectedDEP
        real(rkind) CritError, USTOKES_bar, VSTOKES_bar
        real(rkind) USTOKES_bar_int, VSTOKES_bar_int
        real(rkind) eSum_tot, eSum_tot_int, eWkReal
        real(rkind) eSum_totA, eSum_totA_int
        real(rkind) eSum_totB, eSum_totB_int
        real(rkind) TotalBarotropicErrorUstokes, TotalBarotropicErrorVstokes
        real(rkind) TotalSumUstokes, TotalSumVstokes
        real(rkind) SumHeight
        real(rkind) eJPress_loc, eZetaCorr_loc, eProd, eUint, eVint
# ifndef FIRST_ORDER_ARDHUIN
        eMinMfact=-3
        eMaxMfact=5
# endif
        DO IP=1,MNP
          DO k=1,Nlevel
            z_r(k)=(z_w_wav(IP,k)+z_w_wav(IP,k-1))/2
          END DO
          DO k=0,Nlevel
            z_w_loc(k)=z_w_wav(IP,k)
          END DO
          DO k=2,Nlevel
            PartialU1(k)=(U_wav(IP,k) - U_wav(IP,k-1))/(z_r(k)-z_r(k-1))
            PartialV1(k)=(V_wav(IP,k) - V_wav(IP,k-1))/(z_r(k)-z_r(k-1))
          END DO
          eDep=z_w_loc(Nlevel)-z_w_loc(0)
          PartialU1(1)=PartialU1(2)
          PartialV1(1)=PartialV1(2)
# ifndef FIRST_ORDER_ARDHUIN
          DO k=2,Nlevel-1
            eF1=(z_r(k)-z_r(k-1))/(z_r(k+1)-z_r(k-1))
            eF2=(z_r(k+1)-z_r(k))/(z_r(k+1)-z_r(k-1))
            eDelta=(z_r(k) - z_r(k+1))*(z_r(k-1) - z_r(k))
            PartialU2(k)=(U_wav(IP,k+1)*eF1 + U_wav(IP,k-1)*eF2 - U_wav(IP,k))/eDelta
            PartialV2(k)=(V_wav(IP,k+1)*eF1 + V_wav(IP,k-1)*eF2 - V_wav(IP,k))/eDelta
          END DO
          PartialU2(1)=PartialU2(2)
          PartialV2(1)=PartialV2(2)
          PartialU2(Nlevel)=PartialU2(Nlevel-1)
          PartialV2(Nlevel)=PartialV2(Nlevel-1)
# endif
          eUSTOKES_loc=0
          eVSTOKES_loc=0
          eJpress_loc=0
          eZetaCorr_loc=0
# ifdef FIRST_ORDER_ARDHUIN
          DO IS=1,MSC
            eMult=SPSIG(IS)*DDIR*DS_INCR(IS)
            eWk=WK(IP,IS)
            kD=MIN(KDMAX, eWk*eDep)
            eWkReal=kD/eDep
            eSinh2kd=DSINH(2*kD)
            eSinhkd=DSINH(kD)
            eSinhkd2=eSinhkd**2
            eSigma=SPSIG(IS)
            eUint=0
            eVint=0
            DO ID=1,MDC
              eLoc=AC2(IP,IS,ID)*eMult
              eScal=COSTH(ID)*PartialU1(Nlevel)+SINTH(ID)*PartialV1(Nlevel)
              eZeta=eWk/eSinhkd + (eWk/eSigma)*eScal
              eZetaCorr_loc=eZetaCorr_loc + eLoc*eZeta
              eJPress=G9*(kD/eSinh2kd)*(1/eDep) * eLoc
              eJPress_loc=eJPress_loc + eJPress
              J_PRESSURE(IP)=J_PRESSURE(IP)+eJPress
              eUint=eUint + eLoc*COSTH(ID)
              eVint=eVint + eLoc*SINTH(ID)
            END DO
            DO k=1,Nlevel
              eFrac=(z_r(k) - z_w_loc(0))/eDep
              eHeight=z_w_loc(k)-z_w_loc(k-1)
              eFracB=eHeight/eDep
              eSinc=SINH(kD*eFracB)/(kD*eFracB)
              eQuot1=eSinc*DCOSH(2*kD*eFrac)/eSinhkd2
              eProd=eSigma*eWkReal*eQuot1
              eUSTOKES_loc(k)=eUSTOKES_loc(k) + eUint*eProd
              eVSTOKES_loc(k)=eVSTOKES_loc(k) + eVint*eProd
            ENDDO
          END DO
# else
          DO IS=1,MSC
            eMult=SPSIG(IS)*DDIR*DS_INCR(IS)
            eWk=WK(IP,IS)
            kD=MIN(KDMAX, eWk*eDep)
            eWkReal=kD/eDep
            eSinh2kd=DSINH(2*kD)
            eSinhkd=DSINH(kD)
            eSinhkd2=eSinhkd**2
            eSigma=SPSIG(IS)
            DO ID=1,MDC
              eLoc=AC2(IP,IS,ID)*eMult
              TheInt=0
              DO k=1,Nlevel
                eHeight=z_w_loc(k)-z_w_loc(k-1)
                zMid=0.5*(z_w_loc(k)+z_w_loc(k-1))
                eFrac=(zMid - z_w_loc(0))/eDep
                eFracB=eHeight/eDep
                eSinc=DSINH(eFracB*kD)/(eFracB*kD)
                eQuot=eWkReal*2*DCOSH(2*kD*eFrac)/eSinh2kd
                eFct=U_wav(IP,k)*COSTH(ID)+V_wav(IP,k)*SINTH(ID)
                TheInt=TheInt+eHeight*eFct*eQuot*eSinc
              END DO
              eOmega=eSigma + TheInt*eWkReal
              DO k=1,Nlevel
                MFACT=eSigma/(eOmega - (U_wav(IP,k)*COSTH(ID)+V_wav(IP,k)*SINTH(ID))*eWkReal)
                MFACT=MAX(MFACT, eMinMfact)
                MFACT=MIN(MFACT, eMaxMfact)
                eFrac=(z_r(k) - z_w_loc(0))/eDep
                eHeight=z_w_loc(k)-z_w_loc(k-1)
                eFracB=eHeight/eDep
                eSinc=SINH(kD*eFracB)/(kD*eFracB)
                eQuot1=eSinc*DCOSH(2*kD*eFrac)/eSinhkd2
                USTOKES1=MFACT*eSigma*COSTH(ID)*eWkReal*eQuot1
                VSTOKES1=MFACT*eSigma*SINTH(ID)*eWkReal*eQuot1
                eQuot2=eSinc*DSINH(2*kD*eFrac)/eSinhkd2
                eQuot3=eSinc*(DSINH(kD*eFrac)/eSinhkd)**2
                eScal=PartialU1(k)*COSTH(ID) + PartialV1(k)*SINTH(ID)
                USTOKES2=0.5*(MFACT**2)*COSTH(ID)*eWkReal*eQuot2*eScal
                USTOKES3=0.5*MFACT*PartialU2(k)*eQuot3
                VSTOKES2=0.5*(MFACT**2)*SINTH(ID)*eWkReal*eQuot2*eScal
                VSTOKES3=0.5*MFACT*PartialV2(k)*eQuot3
                USTOKESpart=eLoc*(USTOKES1+USTOKES2+USTOKES3)
                VSTOKESpart=eLoc*(VSTOKES1+VSTOKES2+VSTOKES3)
                eUSTOKES_loc(k)=eUSTOKES_loc(k) + USTOKESpart
                eVSTOKES_loc(k)=eVSTOKES_loc(k) + VSTOKESpart
              ENDDO
              eScal=COSTH(ID)*PartialU1(Nlevel)+SINTH(ID)*PartialV1(Nlevel)
              eZeta=eWk/eSinhkd + (MFACT*eWk/eSigma)*eScal
              eZetaCorr_loc=eZetaCorr_loc + MFACT*eLoc*eZeta
              eJPress=G9*(kD/eSinh2kd)*(1/eDep) * eLoc
              eJPress_loc=eJPress_loc + eJPress
              J_PRESSURE(IP)=J_PRESSURE(IP)+eJPress
            END DO
          END DO
# endif
          DO k=1,Nlevel
            USTOKES_wav(IP,k)=eUSTOKES_loc(k)
            VSTOKES_wav(IP,k)=eVSTOKES_loc(k)
          END DO
          ZETA_CORR(IP)=eZetaCorr_loc
          J_PRESSURE(IP)=eJPress_loc
        ENDDO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE STOKES_STRESS_INTEGRAL_SHYFEM
        USE mod_coupler
        USE DATAPOOL
        implicit none
        integer IP, k, ID, IS
        real(rkind) eF1, eF2, eDelta, TheInt, eDep, eHeight
        real(rkind) eFrac, eFracB, eQuot, TheIntChk
        real(rkind) eFct, eQuot1, eQuot2, eQuot3, eScal, eZeta
        real(rkind) eOmega, eMult, MFACT, kD, eSinc
        real(rkind) USTOKES1, USTOKES2, USTOKES3
        real(rkind) VSTOKES1, VSTOKES2, VSTOKES3
        real(rkind) USTOKESpart, VSTOKESpart, eJPress
        real(rkind) ACLOC, eWk, eSigma, eLoc, eSinhkd, eSinh2kd, eSinhkd2
        real(rkind) zMid, HS, ETOT, MinVal_MFACT, MaxVal_MFACT, MaxVal_eQuot1
        real(rkind) MinVal_MFACT_gl, MaxVal_MFACT_gl, MaxHS, SumHS, AvgHS
        real(rkind) WLM, KLM, AvgStokesNormA, AvgStokesNormB
        real(rkind) PPTAIL, CETAIL, CKTAIL
        real(rkind) ETOT1, EKTOT
        real(rkind) eQuotDispersion, eMaxAC, TotSumAC, eQuotAC, eQuotK
        real(rkind) StokesNormA, StokesNormB, cPhase
        integer IDsel, ISsel, SelectedK
        logical DoTail
        real(rkind) SumNormStokesA(Nlevel), SumNormStokesB(Nlevel)
        real(rkind) SumZetaCorr, MaxZetaCorr, AvgZetaCorr
        real(rkind) eMinMfact, eMaxMfact, SelectedHS
        real(rkind) MaxStokesNorm, MaxValSinc, StokesNorm, SelectedDEP
        real(rkind) CritError, USTOKES_bar, VSTOKES_bar
        real(rkind) USTOKES_bar_int, VSTOKES_bar_int
        real(rkind) eSum_tot, eSum_tot_int, eWkReal
        real(rkind) eSum_totA, eSum_totA_int
        real(rkind) eSum_totB, eSum_totB_int
        real(rkind) TotalBarotropicErrorUstokes, TotalBarotropicErrorVstokes
        real(rkind) TotalSumUstokes, TotalSumVstokes
        real(rkind) SumHeight
        real(rkind) eJPress_loc, eZetaCorr_loc, eProd, eUint, eVint

        DO IP=1,MNP
          eDep=SHYFZETA(NLEV(IP),MNP)
          eUSTOKES_loc=0
          eVSTOKES_loc=0
          eJpress_loc=0
          eZetaCorr_loc=0
          DO IS=1,MSC
            eMult=SPSIG(IS)*DDIR*DS_INCR(IS)
            eWk=WK(IP,IS)
            kD=MIN(KDMAX, eWk*eDep)
            eWkReal=kD/eDep
            eSinh2kd=DSINH(2*kD)
            eSinhkd=DSINH(kD)
            eSinhkd2=eSinhkd**2
            eSigma=SPSIG(IS)
            eUint=0
            eVint=0
            DO ID=1,MDC
              eLoc=AC2(IP,IS,ID)*eMult
              eScal=COSTH(ID)*PartialU1(Nlevel)+SINTH(ID)*PartialV1(Nlevel)
              eZeta=eWk/eSinhkd + (eWk/eSigma)*eScal
              eZetaCorr_loc=eZetaCorr_loc + eLoc*eZeta
              eJPress=G9*(kD/eSinh2kd)*(1/eDep) * eLoc
              eJPress_loc=eJPress_loc + eJPress
              J_PRESSURE(IP)=J_PRESSURE(IP)+eJPress
              eUint=eUint + eLoc*COSTH(ID)
              eVint=eVint + eLoc*SINTH(ID)
            END DO
            DO k=1,Nlevel
              eFrac=(z_r(k) - z_w_loc(0))/eDep
              eHeight=z_w_loc(k)-z_w_loc(k-1)
              eFracB=eHeight/eDep
              eSinc=SINH(kD*eFracB)/(kD*eFracB)
              eQuot1=eSinc*DCOSH(2*kD*eFrac)/eSinhkd2
              eProd=eSigma*eWkReal*eQuot1
              eUSTOKES_loc(k)=eUSTOKES_loc(k) + eUint*eProd
              eVSTOKES_loc(k)=eVSTOKES_loc(k) + eVint*eProd
            ENDDO
          END DO
          DO k=1,Nlevel
            USTOKES_wav(IP,k)=eUSTOKES_loc(k)
            VSTOKES_wav(IP,k)=eVSTOKES_loc(k)
          END DO
          ZETA_CORR(IP)=eZetaCorr_loc
          J_PRESSURE(IP)=eJPress_loc
        ENDDO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PGMCL_ROMS_IN(K,IFILE,IT)
        USE pgmcl_library
        USE datapool
        USE mod_coupler
        implicit none
        INTEGER, INTENT(IN) :: K,IFILE,IT
        integer IP, kLev, i, idx
        real(rkind) u1, v1, u2, v2, z1
# ifdef DEBUG_WWM
        real(rkind) :: MaxUwind, SumUwind, avgUwind
        real(rkind) :: MaxVwind, SumVwind, avgVwind
        real(rkind) :: NbPoint
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: Begin PGMCL_ROMS_IN'
# endif
# ifdef DUMMY_COUPLING
        CALL MPI_INTERP_Dummy_Recv(ArrLocal, 201, OCNid)
# else
        CALL MPI_INTERP_RECV_3D_r8(TheArr_OCNtoWAV_rho, 201, 3, A_wav_rho_3D)
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, After Data receive'
# endif
# ifdef DEBUG_WWM
        MaxUwind=0.0_r8
        SumUwind=0.0_r8
        MaxVwind=0.0_r8
        SumVwind=0.0_r8
        NbPoint=0
# endif
        WATLEVOLD=WATLEV
        DELTAT_WATLEV = MAIN%DTCOUP
        DO idx=1,MNP
          u1=A_wav_rho_3D(idx,1)
          v1=A_wav_rho_3D(idx,2)
# ifdef DEBUG_WWM
          IF (abs(u1).gt.MaxUwind) THEN
            MaxUwind=abs(u1)
          ENDIF
          IF (abs(v1).gt.MaxVwind) THEN
            MaxVwind=abs(v1)
          ENDIF
          SumUwind=SumUwind + abs(u1)
          SumVwind=SumVwind + abs(v1)
          NbPoint=NbPoint+1
# endif
          IP=ReindexPerm_wav(idx)
          u2=u1*DCOS(ang_rho(IP))-v1*DSIN(ang_rho(IP))
          v2=v1*DCOS(ang_rho(IP))+u1*DSIN(ang_rho(IP))
          z1=A_wav_rho_3D(idx,3)
          WINDXY(IP,1)=u2
          WINDXY(IP,2)=v2
          WATLEV(IP)=z1
        END DO
# ifdef DEBUG_WWM
        avgUwind=SumUwind/NbPoint
        avgVwind=SumVwind/NbPoint
        WRITE(DBG%FHNDL,*) 'WAV, MaxUwind=', MaxUwind, ' avgUwind=', avgUwind
        WRITE(DBG%FHNDL,*) 'WAV, MaxVwind=', MaxVwind, ' avgVwind=', avgVwind
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 2'
# endif
# ifdef DUMMY_COUPLING
        CALL MPI_INTERP_Dummy_Recv(ArrLocal, 203, OCNid)
# else
        CALL MPI_INTERP_RECV_3D_r8(TheArr_OCNtoWAV_rho, 203, Nlevel+1, A_wav_rho_3D)
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 3'
# endif
        DO kLev=0,Nlevel
          DO idx=1,MNP
            IP=ReindexPerm_wav(idx)
            z_w_wav(IP,kLev)=A_wav_rho_3D(idx,kLev+1)
          END DO
        END DO
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 4'
# endif
# ifdef DUMMY_COUPLING
        CALL MPI_INTERP_Dummy_Recv(ArrLocal, 204, OCNid)
# else
        CALL MPI_INTERP_RECV_3D_r8(TheArr_OCNtoWAV_u, 204, Nlevel, A_wav_rho_3D)
# endif
        DO kLev=1,Nlevel
          DO idx=1,MNP
            IP=ReindexPerm_wav(idx)
            U_wav(IP,kLev)=A_wav_rho_3D(idx,kLev)
          END DO
        END DO
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 5'
# endif
# ifdef DUMMY_COUPLING
        CALL MPI_INTERP_Dummy_Recv(ArrLocal, 205, OCNid)
# else
        CALL MPI_INTERP_RECV_3D_r8(TheArr_OCNtoWAV_v, 205, Nlevel, A_wav_rho_3D)
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'After the receive'
# endif
        DO kLev=1,Nlevel
          DO idx=1,MNP
            IP=ReindexPerm_wav(idx)
            V_wav(IP,kLev)=A_wav_rho_3D(idx,kLev)
          END DO
        END DO
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 6'
# endif
        DO IP=1,MNP
          DO kLev=1,Nlevel
            u1=U_wav(IP,kLev)
            v1=V_wav(IP,kLev)
            u2=u1*DCOS(ang_rho(IP))-v1*DSIN(ang_rho(IP))
            v2=v1*DCOS(ang_rho(IP))+u1*DSIN(ang_rho(IP))
            U_wav(IP,kLev)=u2
            V_wav(IP,kLev)=v2
          END DO
          CURTXY(IP,1)=u2
          CURTXY(IP,2)=v2
        END DO
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 7'
# endif
      END SUBROUTINE PGMCL_ROMS_IN
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PGMCL_ROMS_OUT(K)
        USE DATAPOOL
        USE pgmcl_library
        USE mod_coupler
# ifdef ST41
        USE W3SRC4MD_OLD, only : Z0, TAUWX, TAUWY, UFRIC, CD
# endif
# ifdef ST42
        USE W3SRC4MD_NEW, only : Z0, TAUWX, TAUWY, UFRIC, CD
# endif
        implicit none
        INTEGER, INTENT(IN)  :: K
        integer IP, kLev, idx
        real(rkind) u1, v1, u2, v2
        real(rkind) HS, TM01, TM02, KLM, WLM, TM10
        real(rkind) UBOT, ORBITAL, BOTEXPER, TMBOT
        real(rkind) FPP, TPP, CPP, WNPP, CGPP, KPP, LPP, PEAKDSPR, PEAKDM, DPEAK
        real(rkind) ETOTS,ETOTC,DM,DSPR
        REAL(RKIND) :: ACLOC(MSC,MDC)
        real(rkind) cPhase, eStokesNorm
        real(rkind) kD
        real(rkind) :: TPPD, KPPD, CGPD, CPPD
# ifdef DEBUG_WWM
        real(rkind) SumNormTau, MaxNormTau, AvgNormTau, eNorm
        real(rkind) AvgUFRICsqr, SumUFRICsqr
        real(rkind) AvgCd, SumCd
        real(rkind) AvgStressCd, SumStressCd, eStressCd, eMag
        real(rkind) AvgAlpha, SumAlpha, eAlpha, NbAlpha
        real(rkind) SumWind, AvgWind
        real(rkind) :: MaxHwave, SumHwave, avgHwave, NbPoint
        real(rkind) :: MaxLwave, SumLwave, avgLwave
        real(rkind) :: MaxStokesNorm, SumStokesNorm, avgStokesNorm
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 1'
# endif
        CALL STOKES_STRESS_INTEGRAL_ROMS
        DO IP=1,MNP
          idx=ReindexPermInv_wav(IP)
          DO kLev=1,Nlevel
            u1=USTOKES_wav(IP,kLev)
            v1=VSTOKES_wav(IP,kLev)
            u2=u1*COS(ang_rho(IP))+v1*SIN(ang_rho(IP))
            v2=v1*COS(ang_rho(IP))-u1*SIN(ang_rho(IP))
            A_wav_rho_3D(idx,kLev)=u2
          END DO
        END DO
# ifdef DUMMY_COUPLING
        CALL MPI_INTERP_Dummy_Send(ArrLocal, 206, OCNid)
        CALL MPI_INTERP_Dummy_Send(ArrLocal, 208, OCNid)
        CALL MPI_INTERP_Dummy_Send(ArrLocal, 210, OCNid)
# else
        CALL MPI_INTERP_SEND_3D_r8(TheArr_WAVtoOCN_rho, 206, Nlevel, A_wav_rho_3D)
        CALL MPI_INTERP_SEND_3D_r8(TheArr_WAVtoOCN_u, 208, Nlevel, A_wav_rho_3D)
        CALL MPI_INTERP_SEND_3D_r8(TheArr_WAVtoOCN_v, 210, Nlevel, A_wav_rho_3D)
# endif
        DO IP=1,MNP
          idx=ReindexPermInv_wav(IP)
          DO kLev=1,Nlevel
            u1=USTOKES_wav(IP,kLev)
            v1=VSTOKES_wav(IP,kLev)
            u2=u1*COS(ang_rho(IP))+v1*SIN(ang_rho(IP))
            v2=v1*COS(ang_rho(IP))-u1*SIN(ang_rho(IP))
            A_wav_rho_3D(idx,kLev)=v2
          END DO
        END DO
# ifdef DUMMY_COUPLING
        CALL MPI_INTERP_Dummy_Send(ArrLocal, 207, OCNid)
        CALL MPI_INTERP_Dummy_Send(ArrLocal, 209, OCNid)
        CALL MPI_INTERP_Dummy_Send(ArrLocal, 211, OCNid)
# else
        CALL MPI_INTERP_SEND_3D_r8(TheArr_WAVtoOCN_rho, 207, Nlevel, A_wav_rho_3D)
        CALL MPI_INTERP_SEND_3D_r8(TheArr_WAVtoOCN_u, 209, Nlevel, A_wav_rho_3D)
        CALL MPI_INTERP_SEND_3D_r8(TheArr_WAVtoOCN_v, 211, Nlevel, A_wav_rho_3D)
# endif
        DO IP=1,MNP
          idx=ReindexPermInv_wav(IP)
          A_wav_rho_3D(idx,1)=J_PRESSURE(IP)
          A_wav_rho_3D(idx,2)=ZETA_CORR(IP)
        END DO
# ifdef DUMMY_COUPLING
        CALL MPI_INTERP_Dummy_Send(ArrLocal, 2047, OCNid)
# else
        CALL MPI_INTERP_SEND_3D_r8(TheArr_WAVtoOCN_rho, 2047, 2, A_wav_rho_3D)
# endif
# ifdef DEBUG_WWM
        SumNormTau=0
        MaxNormTau=0
# endif
        DO IP=1,MNP
          idx=ReindexPermInv_wav(IP)
          u1=TAUWX(IP)
          v1=TAUWY(IP)
          u2=u1*COS(ang_rho(IP))+v1*SIN(ang_rho(IP))
          v2=v1*COS(ang_rho(IP))-u1*SIN(ang_rho(IP))
          A_wav_rho_3D(idx,1)=u2
          A_wav_rho_3D(idx,2)=v2
# ifdef DEBUG_WWM
          eNorm=SQRT(u2*u2 + v2*v2)
          IF (eNorm.gt.MaxNormTau) THEN
            MaxNormTau=eNorm
          END IF
          SumNormTau=SumNormTau + eNorm
# endif
        END DO
# ifdef DEBUG_WWM
        AvgNormTau=SumNormTau / MNP
        WRITE(DBG%FHNDL,*) 'AvgNormTau=', AvgNormTau, 'MaxNormTau=', MaxNormTau
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 5.1'
# endif
# ifdef DUMMY_COUPLING
        CALL MPI_INTERP_Dummy_Send(ArrLocal, 2041, OCNid)
# else
        CALL MPI_INTERP_SEND_3D_r8(TheArr_WAVtoOCN_u, 2041, 2, A_wav_rho_3D)
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 5.2'
# endif
# ifdef DUMMY_COUPLING
        CALL MPI_INTERP_Dummy_Send(ArrLocal, 2039, OCNid)
# else
        CALL MPI_INTERP_SEND_3D_r8(TheArr_WAVtoOCN_v, 2039, 2, A_wav_rho_3D)
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 5.3'
# endif
# ifdef DEBUG_WWM
        SumUFRICsqr=0
        SumCd=0
        SumWind=0
        SumStressCd=0
        SumAlpha=0
        NbAlpha=0
# endif
        DO IP=1,MNP
          idx=ReindexPermInv_wav(IP)
          A_wav_rho_3D(idx,1)=UFRIC(IP)
          A_wav_rho_3D(idx,2)=Z0(IP)
          A_wav_rho_3D(idx,3)=CD(IP)
# ifdef DEBUG_WWM
          SumUFRICsqr=SumUFRICsqr + UFRIC(IP)*UFRIC(IP)
          eMag=SQRT(WINDXY(IP,1)**2 + WINDXY(IP,2)**2)
          eStressCd=CD(IP)*eMag*eMag
          SumWind=SumWind + eMag
          SumStressCd=SumStressCd + eStressCd
          IF (UFRIC(IP).gt.0) THEN
            eAlpha=G9*Z0(IP)/(UFRIC(IP) * UFRIC(IP))
            SumAlpha=SumAlpha + eAlpha
            NbAlpha=NbAlpha+1
          END IF
# endif
        END DO
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 5.4'
# endif
# ifdef DEBUG_WWM
        AvgUFRICsqr=SumUFRICsqr/MNP
        AvgStressCd=SumStressCd/MNP
        AvgAlpha=SumAlpha/NbAlpha
        AvgWind=SumWind/MNP
        AvgCd=SumCd/MNP
        WRITE(DBG%FHNDL,*) 'AvgNormTau=', AvgNormTau, 'AvgNormFV2=', AvgUFRICsqr
        WRITE(DBG%FHNDL,*) 'AvgNormTau=', AvgNormTau, 'AvgCdU2=', AvgStressCd
        WRITE(DBG%FHNDL,*) 'AvgCd=', AvgCd, ' AvgAlpha=', AvgAlpha
        WRITE(DBG%FHNDL,*) 'AvgWind=', AvgWind
# endif
# ifdef DUMMY_COUPLING
        CALL MPI_INTERP_Dummy_Send(ArrLocal, 1453, OCNid)
# else
        CALL MPI_INTERP_SEND_3D_r8(TheArr_WAVtoOCN_rho, 1453, 3, A_wav_rho_3D)
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 5.5'
# endif
# ifdef DEBUG_WWM
        MaxHwave=0.0
        SumHwave=0.0
        MaxLwave=0.0
        SumLwave=0.0
        SumStokesNorm=0
        MaxStokesNorm=0
        NbPoint=0.0
# endif
        DO IP = 1, MNP
          idx=ReindexPermInv_wav(IP)
          ACLOC = AC2(IP,:,:)
          CALL MEAN_PARAMETER(IP,ACLOC,MSC,HS,TM01,TM02,TM10,KLM,WLM)
# ifdef DEBUG_WWM
          IF (HS.gt.MaxHwave) THEN
            MaxHwave=HS
          ENDIF
          SumHwave=SumHwave + HS
          IF (WLM.gt.MaxLwave) THEN
            MaxLwave=WLM
          ENDIF
          SumLwave=SumLwave + WLM
          kD=MIN(KDMAX, KLM*DEP(IP))
          cPhase=SQRT((G9/KLM)*REAL(DSINH(kD)/DCOSH(kD)) )
          eStokesNorm=(G9*HS*HS/REAL(16))*2*(KLM/cPhase)
          IF (eStokesNorm.ne.eStokesNorm) THEN
            WRITE(DBG%FHNDL,*) 'eStokesNorm=', eStokesNorm
            WRITE(DBG%FHNDL,*) 'KLM=', KLM, 'WLM=', WLM
            WRITE(DBG%FHNDL,*) 'cPhase=', cPhase, 'kD=', kD
            WRITE(DBG%FHNDL,*) 'HS=', HS, ' DEP=', DEP(IP)
          END IF
          IF (eStokesNorm.gt.MaxStokesNorm) THEN
            MaxStokesNorm=eStokesNorm
          END IF
          SumStokesNorm=SumStokesNorm + eStokesNorm
          NbPoint=NbPoint + 1
# endif
          A_wav_rho_3D(idx,1)=HS
          A_wav_rho_3D(idx,2)=TM01
          A_wav_rho_3D(idx,3)=TM02
          A_wav_rho_3D(idx,4)=KLM
          A_wav_rho_3D(idx,5)=WLM
        END DO
# ifdef DEBUG_WWM
        avgHwave=SumHwave/NbPoint
        avgLwave=SumLwave/NbPoint
        avgStokesNorm=SumStokesNorm/NbPoint
        WRITE(DBG%FHNDL,*) 'WAV, MaxHwave=', MaxHwave, ' avgHwave=', avgHwave
        WRITE(DBG%FHNDL,*) 'WAV, MaxLwave=', MaxLwave, ' avgLwave=', avgLwave
        WRITE(DBG%FHNDL,*) 'WAV, MaxStokesNorm=', MaxStokesNorm, ' avgStokesNorm=', avgStokesNorm
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 6'
# endif
# ifdef DUMMY_COUPLING
        CALL MPI_INTERP_Dummy_Send(ArrLocal, 212, OCNid)
# else
        CALL MPI_INTERP_SEND_3D_r8(TheArr_WAVtoOCN_rho, 212, 5, A_wav_rho_3D)
# endif

        DO IP = 1, MNP
          idx=ReindexPermInv_wav(IP)
          ACLOC = AC2(IP,:,:)
          CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT)
          A_wav_rho_3D(idx,1)=ORBITAL
          A_wav_rho_3D(idx,2)=TMBOT
          A_wav_rho_3D(idx,3)=DISSIPATION(IP)
          A_wav_rho_3D(idx,4)=QBLOCAL(IP)
        END DO
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 8'
# endif
# ifdef DUMMY_COUPLING
        CALL MPI_INTERP_Dummy_Send(ArrLocal, 213, OCNid)
# else
        CALL MPI_INTERP_SEND_3D_r8(TheArr_WAVtoOCN_rho, 213, 4, A_wav_rho_3D)
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 9'
# endif
        DO IP = 1, MNP
          idx=ReindexPermInv_wav(IP)
          ACLOC = AC2(IP,:,:)
          CALL MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,MSC,ETOTS,ETOTC,DM,DSPR)
          CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
          A_wav_rho_3D(idx,1)=DM
          A_wav_rho_3D(idx,2)=TPP
          A_wav_rho_3D(idx,3)=DSPR
          A_wav_rho_3D(idx,4)=PEAKDSPR
          A_wav_rho_3D(idx,5)=PEAKDM
        END DO
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 10'
# endif
# ifdef DUMMY_COUPLING
        CALL MPI_INTERP_Dummy_Send(ArrLocal, 214, OCNid)
# else
        CALL MPI_INTERP_SEND_3D_r8(TheArr_WAVtoOCN_rho, 214, 5, A_wav_rho_3D)
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 11'
# endif
      END SUBROUTINE
#endif
