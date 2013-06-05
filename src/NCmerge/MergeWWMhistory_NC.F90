      PROGRAM MERGE_WWM_HISTORY_NC
      USE NETCDF
      USE DATAPOOL
      IMPLICIT NONE
      CHARACTER(LEN = 100) :: HisPrefix
      CHARACTER(LEN = 100) :: VarName
      CHARACTER(LEN = 100) :: BinFile
      character(len = 1000) :: FILE_NAME_SPLIT
      character(len = 1000) :: FILE_NAME_MERGE
      character(len = 728) :: CHRERR
      character(len=40) :: eStr, eStrUnit
      character(len=80) :: eStrFullName
      character (len = *), parameter :: UNITS = "units"
      character (len = *), parameter :: FULLNAME = "full-name"
      integer, dimension(nf90_max_var_dims) :: dimids

      integer nbArg, nbProc, iProc, FHNDL
      integer sizBlock
      integer IP, len, nb1, nb2, eMNP
      integer iret, idim_id, ncid, ncidB, var_id, itime_id
      integer nnode_dims, ntime_dims
      integer iTime, nbTime, iEnt
      integer len1, len2
      character (len = *), parameter :: CallFct="MERGE_WWM_NC_FILE"

      integer, allocatable :: ListMNP(:)
      integer, allocatable :: iplg(:)
      integer ifile, I
      integer :: eInt(1)
     
      logical test
      real*4, allocatable :: TheOut_r4(:,:)
      real*8, allocatable :: TheOut_r8(:,:)
      real*4, allocatable :: TheR_r4(:)
      real*8, allocatable :: TheR_r8(:)
      real(rkind), allocatable :: LTimeDay(:)
      real(rkind) :: eTimeDay

      nbArg=command_argument_count()
      IF (nbArg .ge. 2) THEN
        PRINT *, 'MergeWWMhistory_NC merges netcdf history files'
        Print *, 'It is supposed to be used just like WWM'
        Print *, ''
        Print *, 'i.e. MergeWWMhistory_NC'
        Print *, 'or   MergeWWMhistory_NC wwminput.nml'
        Print *, 'with the wwminput the same as the one of the run'
        STOP
      END IF
      IF (nbArg.eq.0) THEN
        INP%FNAME  = 'wwminput.nml'
      ELSE
        CALL GET_COMMAND_ARGUMENT(1, INP%FNAME)
      ENDIF
!
! Determine nbTime and np_global
! Also read ocean_time_day
!

!      INP%FHNDL=440
!      open(INP%FHNDL, FILE=TRIM(INP%FNAME))
!      GRD%FHNDL=230
!      STAT%FHNDL=235
!      open(STAT%FHNDL, FILE='wwmstat.out', status='unknown')
      CALL INIT_FILE_HANDLES()
      CALL READ_WWMINPUT
!      CALL INIT_ARRAYS
      CLOSE(INP%FHNDL)
      open(GRD%FHNDL, FILE=TRIM(GRD%FNAME))
!
      Print *, 'GRD%FNAME=', TRIM(GRD%FNAME)

      ALLOCATE(XP(MNP))
      ALLOCATE(YP(MNP))
      ALLOCATE(DEP(MNP))
      ALLOCATE(INVSPHTRANS(MNP,2))
      ALLOCATE(INE(3,MNE))
      ALLOCATE(IEN(6,MNE))
      ALLOCATE(TRIA(MNE))
!
      CALL READ_SPATIAL_GRID
      CALL SPATIAL_GRID
      close(GRD%FHNDL)
!
! Creating the first FILE_NAME
!
      ifile=1
      iProc=1
      CALL GET_FILE_NAME_SPLIT(FILE_NAME_SPLIT, ifile, iProc)
!
! Reading nproc
!
      Print *, 'FILE_NAME_SPLIT=', TRIM(FILE_NAME_SPLIT)
      iret = NF90_OPEN(TRIM(FILE_NAME_SPLIT), NF90_NOWRITE, ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, iret)

      iret = nf90_inq_varid(ncid, 'nproc', var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, iret)

      iret = NF90_GET_VAR(ncid, var_id, eInt)
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, iret)
      nproc=eInt(1)

      iret = NF90_CLOSE(ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 4, iret)
!
! Read the MNP
!
      allocate(ListMNP(nproc))
      DO iProc=1,nproc
        ifile=1
        CALL GET_FILE_NAME_SPLIT(FILE_NAME_SPLIT, ifile, iProc)
        iret = NF90_OPEN(TRIM(FILE_NAME_SPLIT), NF90_NOWRITE, ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 5, iret)
        iret = nf90_inq_dimid(ncid, 'mnp', idim_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 6, iret)
        iret = nf90_inquire_dimension(ncid, idim_id, len = eMNP)
        CALL GENERIC_NETCDF_ERROR(CallFct, 7, iret)
        ListMNP(iProc)=eMNP
        Print *, 'iProc=', iProc, ' eMNP=', eMNP
        iret = NF90_CLOSE(ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 8, iret)
      END DO
!
! Now looping ...
!
      ifile=1
      IF (NF90_RUNTYPE == NF90_REAL) THEN
        allocate(TheOut_r4(OUTVARS_COMPLETE, MNP))
      ELSE
        allocate(TheOut_r8(OUTVARS_COMPLETE, MNP))
      END IF
      DO
        iProc=1
        CALL GET_FILE_NAME_SPLIT(FILE_NAME_SPLIT, ifile, iProc)
        INQUIRE(FILE=TRIM(FILE_NAME_SPLIT),EXIST=test)
        IF (test .eqv. .FALSE.) THEN
          EXIT
        END IF
        Print *, 'ifile=', ifile
        !
        ! Reading the time
        !
        iret = nf90_open(TRIM(FILE_NAME_SPLIT), NF90_NOWRITE, ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 9, iret)
        iret = nf90_inq_varid(ncid, 'ocean_time_day', itime_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 10, iret)
        iret = NF90_INQUIRE_VARIABLE(ncid, ITIME_ID, dimids = dimids)
        CALL GENERIC_NETCDF_ERROR(CallFct, 11, iret)
        iret = nf90_inquire_dimension(ncid, dimids(1), len = nbTime)
        CALL GENERIC_NETCDF_ERROR(CallFct, 12, iret)
        ALLOCATE(LTimeDay(nbTime))
        iret = NF90_GET_VAR(ncid, itime_id, LTimeDay)
        CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
        iret=nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 14, iret)
        !
        ! Now creating the merged file.
        !
        CALL GET_FILE_NAME_MERGE(FILE_NAME_MERGE, ifile)
        Print *, 'FILE_NAME_MERGE=', TRIM(FILE_NAME_MERGE)
        iret = nf90_create(TRIM(FILE_NAME_MERGE), NF90_CLOBBER, ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 15, iret)
        CALL WRITE_NETCDF_HEADERS_1(ncid, -1, MULTIPLEOUT_HIS, MNP, MNE)
        iret=nf90_inq_dimid(ncid, 'mnp', nnode_dims)
        CALL GENERIC_NETCDF_ERROR(CallFct, 16, iret)
        iret=nf90_inq_dimid(ncid, 'ocean_time', ntime_dims)
        CALL GENERIC_NETCDF_ERROR(CallFct, 17, iret)
        DO I=1,OUTVARS_COMPLETE
          IF (VAROUT_HISTORY%LVAR(I)) THEN
            CALL NAMEVARIABLE(I, eStr, eStrFullName, eStrUnit)
            iret=nf90_def_var(ncid,TRIM(eStr),NF90_OUTTYPE_HIS,(/ nnode_dims, ntime_dims /),var_id)
            CALL GENERIC_NETCDF_ERROR(CallFct, 18, iret)
            iret=nf90_put_att(ncid,var_id,UNITS,TRIM(eStrUnit))
            CALL GENERIC_NETCDF_ERROR(CallFct, 19, iret)
            iret=nf90_put_att(ncid,var_id,FULLNAME,TRIM(eStrFullName))
            CALL GENERIC_NETCDF_ERROR(CallFct, 20, iret)
          END IF
        END DO
        iret=nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 21, iret)
        !
        iret=nf90_open(TRIM(FILE_NAME_MERGE), NF90_WRITE, ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 22, iret)
        WriteOutputProcess_his=.TRUE.
        CALL WRITE_NETCDF_HEADERS_2(ncid, MULTIPLEOUT_HIS,      &
    &       WriteOutputProcess_his, MNP, MNE)
        iret=nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 23, iret)
        !
        ! Now reading the number of times we need to work
        !
        iret = nf90_open(TRIM(FILE_NAME_MERGE), NF90_WRITE, ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 24, iret)
        iret = nf90_inq_varid(ncid, 'ocean_time_day', itime_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 25, iret)
        DO iTime=1,nbTime
!          Print *, '  iTime=', iTime, '/', nbTime
          eTimeDay=LTimeDay(iTime)
          CALL WRITE_NETCDF_TIME(ncid, iTime, eTimeDay)
          DO iProc=1,nproc
            eMNP=ListMNP(iProc)
            IF (NF90_RUNTYPE == NF90_REAL) THEN
              allocate(TheR_r4(eMNP))
            ELSE
              allocate(TheR_r8(eMNP))
            END IF
            ALLOCATE(iplg(eMNP))


            CALL GET_FILE_NAME_SPLIT(FILE_NAME_SPLIT, ifile, iProc)
            CALL GENERIC_NETCDF_ERROR(CallFct, 26, iret)
            iret = nf90_open(TRIM(FILE_NAME_SPLIT), NF90_NOWRITE, ncidB)
            CALL GENERIC_NETCDF_ERROR(CallFct, 27, iret)
            iret = nf90_inq_varid(ncidB, 'iplg', var_id)
            CALL GENERIC_NETCDF_ERROR(CallFct, 28, iret)
            iret = NF90_GET_VAR(ncidB, var_id, iplg)
            CALL GENERIC_NETCDF_ERROR(CallFct, 29, iret)
            DO I=1,OUTVARS_COMPLETE
              IF (VAROUT_HISTORY%LVAR(I)) THEN
                CALL NAMEVARIABLE(I, eStr, eStrFullName, eStrUnit)
                iret = nf90_inq_varid(ncidB, TRIM(eStr), var_id)
                CALL GENERIC_NETCDF_ERROR(CallFct, 30, iret)
                IF (NF90_RUNTYPE == NF90_REAL) THEN
                  iret = NF90_GET_VAR(ncidB, var_id, TheR_r4, start=(/1, iTime/), count=(/eMNP, 1/))
                ELSE
                  iret = NF90_GET_VAR(ncidB, var_id, TheR_r8, start=(/1, iTime/), count=(/eMNP, 1/))
                END IF
                CALL GENERIC_NETCDF_ERROR(CallFct, 31, iret)
                IF (NF90_RUNTYPE == NF90_REAL) THEN
                  DO IP=1,eMNP
                    TheOut_r4(I,iplg(IP))=TheR_r4(IP)
                  END DO
                ELSE
                  DO IP=1,eMNP
                    TheOut_r8(I,iplg(IP))=TheR_r8(IP)
                  END DO
                END IF
              END IF
            END DO
            iret=nf90_close(ncidB)
            CALL GENERIC_NETCDF_ERROR(CallFct, 32, iret)
            IF (NF90_RUNTYPE == NF90_REAL) THEN
              deallocate(TheR_r4)
            ELSE
              deallocate(TheR_r8)
            END IF
            deallocate(iplg)
          END DO
          DO I=1,OUTVARS_COMPLETE
            IF (VAROUT_HISTORY%LVAR(I)) THEN
              CALL NAMEVARIABLE(I, eStr, eStrFullName, eStrUnit)
              iret = nf90_inq_varid(ncid, TRIM(eStr), var_id)
              CALL GENERIC_NETCDF_ERROR(CallFct, 33, iret)
              IF (NF90_RUNTYPE == NF90_REAL) THEN
                iret = NF90_PUT_VAR(ncid, var_id, TheOut_r4(I,:), start=(/1, iTime/), count=(/MNP, 1/))
              ELSE
                iret = NF90_PUT_VAR(ncid, var_id, TheOut_r8(I,:), start=(/1, iTime/), count=(/MNP, 1/))
              ENDIF
              CALL GENERIC_NETCDF_ERROR(CallFct, 34, iret)
            END IF
          END DO
        END DO
        iret=nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 35, iret)
        DEALLOCATE(LTimeDay)
        ifile=ifile+1
      END DO
  10  FORMAT (a,i4.4,'.nc')
      END PROGRAM
!**********************************************************************
!*                                                                    *
!**********************************************************************    
      SUBROUTINE GET_FILE_NAME_SPLIT(FILE_NAME, ifile, iProc)
      USE DATAPOOL
      implicit none
      integer, intent(in) :: ifile, iProc
      character(len=1000), intent(out) :: FILE_NAME
      character(len =256) :: PRE_FILE_NAME
      integer LPOS, POSITION_BEFORE_POINT
      LPOS=POSITION_BEFORE_POINT(OUT_HISTORY%FNAME)
      IF (OUT_HISTORY%IDEF.gt.0) THEN
         WRITE (PRE_FILE_NAME,10) OUT_HISTORY%FNAME(1:LPOS),ifile
  10     FORMAT (a,'_',i4.4)
      ELSE
         WRITE (PRE_FILE_NAME,20) OUT_HISTORY%FNAME(1:LPOS)
  20     FORMAT (a)
      ENDIF
      WRITE (FILE_NAME,40) TRIM(PRE_FILE_NAME),iProc-1
  40  FORMAT (a,'_',i4.4,'.nc')
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************    
      SUBROUTINE GET_FILE_NAME_MERGE(FILE_NAME, ifile)
      USE DATAPOOL
      implicit none
      integer, intent(in) :: ifile
      character(len=1000), intent(out) :: FILE_NAME
      integer LPOS, POSITION_BEFORE_POINT
      LPOS=POSITION_BEFORE_POINT(OUT_HISTORY%FNAME)
      IF (OUT_HISTORY%IDEF.gt.0) THEN
         WRITE (FILE_NAME,10) OUT_HISTORY%FNAME(1:LPOS),ifile
  10     FORMAT (a,'_',i4.4,'.nc')
      ELSE
         WRITE (FILE_NAME,20) OUT_HISTORY%FNAME(1:LPOS)
  20     FORMAT (a,'nc')
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_FILE_HANDLES()
#ifdef MPI_PARALL_GRID
         use elfe_msgp
#endif
         USE DATAPOOL
         IMPLICIT NONE
         CHARACTER (LEN = 30) :: FDB
         INTEGER              :: LFDB

#ifdef SELFE
            INP%FNAME  = 'wwminput.nml'
#endif
            CHK%FNAME  = 'wwmcheck.nml'
           QSTEA%FNAME = 'qstea.out'
         WINDBG%FNAME  = 'winddbg.out'
         IOBPOUT%FNAME = 'iobp.out'
        IOBPDOUT%FNAME = 'iobpd.out'
!
!2do ... dinstinguish between binary and ascii stuff ...
!
             BND%FHNDL  = STARTHNDL + 1
             WIN%FHNDL  = STARTHNDL + 2
             CUR%FHNDL  = STARTHNDL + 3
             WAT%FHNDL  = STARTHNDL + 4
             WAV%FHNDL  = STARTHNDL + 5
             CHK%FHNDL  = STARTHNDL + 6
           HOTIN%FHNDL  = STARTHNDL + 7
          HOTOUT%FHNDL  = STARTHNDL + 8
             INP%FHNDL  = STARTHNDL + 9
             GRD%FHNDL  = STARTHNDL + 10
          GRDCOR%FHNDL  = STARTHNDL + 11

           QSTEA%FHNDL  = STARTHNDL + 12

         IOBPOUT%FHNDL  = STARTHNDL + 13
         IOBPDOUT%FHNDL = STARTHNDL + 14

         DBG%FHNDL      = STARTHNDL + 15 
         STAT%FHNDL     = STARTHNDL + 16 
         WINDBG%FHNDL   = STARTHNDL + 17 

#ifndef MPI_PARALL_GRID
         open(DBG%FHNDL,file='wwmdbg.out',status='unknown') !non-fatal errors
         open(STAT%FHNDL,file='wwmstat.out',status='unknown') !non-fatal errors
         open(WINDBG%FHNDL,file='windbg.out',status='unknown') !non-fatal errors
#else
# ifdef SELFE
         FDB  ='wwmdbg_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(DBG%FHNDL,file='outputs/'//fdb,status='replace') 
         FDB  ='wwmstat_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(STAT%FHNDL,file='outputs/'//fdb,status='replace') 
         FDB  ='windbg_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(WINDBG%FHNDL,file='outputs/'//fdb,status='replace') 
# else
         FDB  ='wwmdbg_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(DBG%FHNDL,file=fdb,status='replace') 
         FDB  ='wwmstat_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(STAT%FHNDL,file=fdb,status='replace') 
         FDB  ='windbg_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(WINDBG%FHNDL,file=fdb,status='replace')
# endif
#endif
         CALL TEST_FILE_EXIST_DIE("Missing input file : ", TRIM(INP%FNAME))
#ifdef MPIP_PARALL_GRID
         IF (myrank == 0) THEN
#endif
         WRITE(STAT%FHNDL,*) 'Input Filename   =', TRIM(INP%FNAME)
         WRITE(STAT%FHNDL,*) 'Check Filename   =', TRIM(CHK%FNAME)
         WRITE(STAT%FHNDL,*) 'Qstea Filename   =', TRIM(QSTEA%FNAME)
         WRITE(STAT%FHNDL,*) 'Iobp Filename    =', TRIM(IOBPOUT%FNAME)
         WRITE(STAT%FHNDL,*) 'Iobpd Filename   =', TRIM(IOBPDOUT%FNAME)
         WRITE(STAT%FHNDL,*) 'WindDbg Filename =', TRIM(WINDBG%FNAME)
#ifdef MPIP_PARALL_GRID
         ENDIF
#endif

         OPEN( INP%FHNDL,      FILE = TRIM(INP%FNAME))
         OPEN( CHK%FHNDL,      FILE = TRIM(CHK%FNAME))
         OPEN( QSTEA%FHNDL,    FILE = TRIM(QSTEA%FNAME))
         OPEN( IOBPOUT%FHNDL,  FILE = TRIM(IOBPOUT%FNAME))
         OPEN( IOBPDOUT%FHNDL, FILE = TRIM(IOBPDOUT%FNAME))

           OUT1D%FHNDL = STARTHNDL + 18 
            MISC%FHNDL = STARTHNDL + 19 
         OUTSP1D%FHNDL = STARTHNDL + 20 
         OUTPARM%FHNDL = STARTHNDL + 21 
         OUTSP2D%FHNDL = STARTHNDL + 22 

         OUT%FHNDL     = STARTHNDL + 23 

      END SUBROUTINE
