      PROGRAM MERGE_WWM_HOTFILE_NC
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
      integer nnode_dims, ntime_dims, msc_dims, mdc_dims, ac_id
      integer iTime, nbTime, iEnt
      integer len1, len2
      character (len = *), parameter :: CallFct="MERGE_WWM_NC_FILE"

      integer, allocatable :: ListMNP(:)
      integer, allocatable :: iplg(:)
      integer I, nproc
      integer :: eInt(1)
     
      logical test
      real(rkind), allocatable :: AC(:,:,:)
      real(rkind), allocatable :: ACloc(:,:,:)
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
      CLOSE(INP%FHNDL)
      open(GRD%FHNDL, FILE=TRIM(GRD%FNAME))
!
      Print *, 'GRD%FNAME=', TRIM(GRD%FNAME)

      ALLOCATE(XP(MNP))
      ALLOCATE(YP(MNP))
      ALLOCATE(DEP(MNP))
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
      iProc=1
      CALL GET_FILE_NAME_SPLIT(FILE_NAME_SPLIT, iProc)
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
        CALL GET_FILE_NAME_SPLIT(FILE_NAME_SPLIT, iProc)
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
      allocate(AC(MNP,MSC,MDC))
      iProc=1
      CALL GET_FILE_NAME_SPLIT(FILE_NAME_SPLIT, iProc)
      INQUIRE(FILE=TRIM(FILE_NAME_SPLIT),EXIST=test)
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
      MULTIPLEOUT_HOT=0
      WriteOutputProcess_hot=.TRUE.
      CALL GET_FILE_NAME_MERGE(FILE_NAME_MERGE)
      Print *, 'FILE_NAME_MERGE=', TRIM(FILE_NAME_MERGE)
      iret = nf90_create(TRIM(FILE_NAME_MERGE), NF90_CLOBBER, ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 15, iret)

      CALL WRITE_NETCDF_HEADERS_1(ncid, nbTime, MULTIPLEOUT_HOT, MNP, MNE)

      iret=nf90_inq_dimid(ncid, 'mnp', nnode_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 16, iret)

      iret=nf90_inq_dimid(ncid, 'ocean_time', ntime_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 17, iret)

      iret=nf90_inq_dimid(ncid, 'msc', msc_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 18, iret)

      iret=nf90_inq_dimid(ncid, 'mdc', mdc_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 19, iret)

      iret=nf90_def_var(ncid,'ac',NF90_RUNTYPE,(/ nnode_dims, msc_dims, mdc_dims, ntime_dims/),ac_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 20, iret)

      iret=nf90_put_att(ncid,ac_id,UNITS,'unknown')
      CALL GENERIC_NETCDF_ERROR(CallFct, 21, iret)

      iret=nf90_close(ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 22, iret)
      !
      ! Putting the header
      !
      iret = nf90_open(TRIM(FILE_NAME_MERGE), nf90_write, ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 23, iret)
      CALL WRITE_NETCDF_HEADERS_2(ncid, MULTIPLEOUT_HOT, WriteOutputProcess_hot, MNP, MNE)
      iret=nf90_close(ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 24, iret)
      !
      ! Now reading data and writing it.
      !
      iret = nf90_open(TRIM(FILE_NAME_MERGE), NF90_WRITE, ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 25, iret)

      iret = nf90_inq_varid(ncid, 'ocean_time_day', itime_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 26, iret)
      DO iTime=1,nbTime
        Print *, '  iTime=', iTime, '/', nbTime
        eTimeDay=LTimeDay(iTime)
        CALL WRITE_NETCDF_TIME(ncid, iTime, eTimeDay)
        DO iProc=1,nproc
          eMNP=ListMNP(iProc)
          allocate(ACloc(eMNP, MSC, MDC))
          ALLOCATE(iplg(eMNP))
          CALL GET_FILE_NAME_SPLIT(FILE_NAME_SPLIT, iProc)
          iret = nf90_open(TRIM(FILE_NAME_SPLIT), NF90_NOWRITE, ncidB)
          CALL GENERIC_NETCDF_ERROR(CallFct, 27, iret)
          iret = nf90_inq_varid(ncidB, 'iplg', var_id)
          CALL GENERIC_NETCDF_ERROR(CallFct, 28, iret)
          iret = NF90_GET_VAR(ncidB, var_id, iplg)
          CALL GENERIC_NETCDF_ERROR(CallFct, 29, iret)

          iret = NF90_GET_VAR(ncidB, var_id, ACloc, start=(/1, 1, 1, iTime/), count=(/eMNP, MSC, MDC, 1/))
          DO IP=1,eMNP
            AC(iplg(IP),:,:)=ACloc(IP,:,:)
          END DO
          iret=nf90_close(ncidB)
          CALL GENERIC_NETCDF_ERROR(CallFct, 30, iret)
          deallocate(iplg)
          deallocate(ACloc)
        END DO
        !
        iret = nf90_inq_varid(ncid, "ac", var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 31, iret)
        iret = NF90_PUT_VAR(ncid, var_id, AC, start=(/1, 1, 1, iTime/), count=(/MNP, MSC, MDC, 1/))
        CALL GENERIC_NETCDF_ERROR(CallFct, 32, iret)
      END DO
      iret=nf90_close(ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 33, iret)
      END PROGRAM
!**********************************************************************
!*                                                                    *
!**********************************************************************    
      SUBROUTINE GET_FILE_NAME_SPLIT(FILE_NAME, iProc)
      USE DATAPOOL
      implicit none
      integer, intent(in) :: iProc
      integer, parameter :: powerproc = 6
      character(len=1000), intent(out) :: FILE_NAME
      character(len=6) :: eStrProc
      integer LPOS, POSITION_BEFORE_POINT
      LPOS=POSITION_BEFORE_POINT(HOTOUT%FNAME)
      CALL GETSTRING(powerproc, iProc, eStrProc)
      WRITE (FILE_NAME,50) TRIM(HOTOUT%FNAME(1:LPOS)),eStrProc
  50  FORMAT (a,'_',a,'.nc')
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************    
      SUBROUTINE GET_FILE_NAME_MERGE(FILE_NAME)
      USE DATAPOOL
      implicit none
      character(len=1000), intent(out) :: FILE_NAME
      integer LPOS, POSITION_BEFORE_POINT
      LPOS=POSITION_BEFORE_POINT(HOTOUT%FNAME)
      WRITE (FILE_NAME,50) TRIM(HOTOUT%FNAME(1:LPOS))
  50  FORMAT (a,'.nc')
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
!**********************************************************************
!*                                                                    *
!**********************************************************************    
