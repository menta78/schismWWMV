      PROGRAM MERGE_WWM_NC_FILE
      USE NETCDF
      IMPLICIT NONE
      CHARACTER(LEN = 1000) :: NcFile
      CHARACTER(LEN = 1000) :: BinFile
      character(len = 728) :: CHRERR
      integer, dimension(nf90_max_var_dims) :: dimids

      integer nbArg, nbProc, iProc, FHNDL
      integer sizBlock
      integer IP, len, nb1, nb2, eMNP
      integer ISTAT, idim_id, ncid, var_id, itime_id
      integer mne_dims, mnp_dims, ntime_dims, one_dims, three_dims
      integer iret, varHs_id, varTime_id
      integer iTime, nbTime, iEnt
      integer np_global, ne_global, idx, len1, len2
      integer i, j, k, iegb
      REAL*8, PARAMETER   :: DAY2SEC  = 86400

      real*4 :: TIME_4
      real*8, allocatable :: XP(:), YP(:), DEP(:)
      real*4, allocatable :: TheMat(:,:)
      real*4, allocatable :: TheR(:,:)
      real*4, allocatable :: TheWrite(:)
      real*8, allocatable :: LTimeDay(:)
      real*4, allocatable :: Aread(:), Bread(:), Cread(:)
      integer, allocatable :: INE(:,:)
      real*8 :: eTime(1)
      real*4 :: avgHs

      nbArg=command_argument_count()
      IF (nbArg == 0) THEN
        PRINT *, 'Reading from a bin file'
        Print *, 'and create a netcdf file'
        Print *, 'ConvertBin_to_NC [BinFile] [NcFile]'
        Print *, ''
        Print *, 'We also need the hgrid.ll file'
        STOP
      END IF
      CALL GET_COMMAND_ARGUMENT(1, BinFile)
      CALL GET_COMMAND_ARGUMENT(2, NcFile)
!
! Reading the grid file
!
      open(14,file='hgrid.ll',status='old')
      read(14,*); read(14,*) ne_global,np_global
      Print *, 'np_global=', np_global, 'ne_global=', ne_global
      allocate(XP(np_global))
      allocate(YP(np_global))
      allocate(DEP(np_global))
      allocate(INE(3,ne_global))
      do i=1,np_global
        read(14,*) idx,XP(i),YP(i), DEP(i)
      enddo
      Print *, 'After coordinate read'
      do i=1,ne_global
        read(14,*) iegb,j,(INE(k,iegb),k=1,3)
        if(j/=3) then
          Print *, 'Error in file reading'
          STOP
        endif
      enddo
      Print *, 'After triangle read'
      close(14)
!
! Determine nbTime and np_global
! Also read ocean_time_day
!
      Print *, 'NcFile=', TRIM(NcFile)
      iret = NF90_CREATE(TRIM(NcFile), NF90_CLOBBER, ncid)
      CALL GENERIC_NETCDF_ERROR(1, iret)

      iret = nf90_def_dim(ncid, 'one', 1, one_dims)
      CALL GENERIC_NETCDF_ERROR(2, iret)

      iret = nf90_def_dim(ncid, 'three', 3, three_dims)
      CALL GENERIC_NETCDF_ERROR(3, iret)

      iret = nf90_def_dim(ncid, 'mnp', np_global, mnp_dims)
      CALL GENERIC_NETCDF_ERROR(4, iret)

      iret = nf90_def_dim(ncid, 'mne', ne_global, mne_dims)
      CALL GENERIC_NETCDF_ERROR(5, iret)

      iret = nf90_def_dim(ncid, 'ocean_time', NF90_UNLIMITED, ntime_dims)
      CALL GENERIC_NETCDF_ERROR(6, iret)

      iret=nf90_def_var(ncid,'lon',NF90_DOUBLE,(/mnp_dims/),var_id)
      CALL GENERIC_NETCDF_ERROR(7, iret)

      iret=nf90_def_var(ncid,'lat',NF90_DOUBLE,(/mnp_dims/),var_id)
      CALL GENERIC_NETCDF_ERROR(8, iret)

      iret=nf90_def_var(ncid,'ele',NF90_INT,(/three_dims, mne_dims/),var_id)
      CALL GENERIC_NETCDF_ERROR(9, iret)
      iret=nf90_def_var(ncid,'ocean_time',NF90_DOUBLE,(/ ntime_dims /),varTime_id)
      CALL GENERIC_NETCDF_ERROR(10, iret)
      iret=nf90_def_var(ncid,'Hwave',NF90_REAL,(/ mnp_dims, ntime_dims /),varHs_id)
      CALL GENERIC_NETCDF_ERROR(11, iret)

      iret = NF90_CLOSE(ncid)
      CALL GENERIC_NETCDF_ERROR(12, iret)
!
! Now writing the data
!
      iret=nf90_open(TRIM(NcFile), NF90_WRITE, ncid)
      CALL GENERIC_NETCDF_ERROR(13, iret)

      iret=nf90_inq_varid(ncid,'lon',var_id)
      CALL GENERIC_NETCDF_ERROR(14, iret)
      iret=nf90_put_var(ncid,var_id,XP,start = (/1/), count=(/np_global/))
      CALL GENERIC_NETCDF_ERROR(15, iret)

      iret=nf90_inq_varid(ncid,'lat',var_id)
      CALL GENERIC_NETCDF_ERROR(16, iret)
      iret=nf90_put_var(ncid,var_id,YP,start = (/1/), count=(/np_global/))
      CALL GENERIC_NETCDF_ERROR(17, iret)


      iret=nf90_inq_varid(ncid,'ele',var_id)
      CALL GENERIC_NETCDF_ERROR(18, iret)
      iret=nf90_put_var(ncid,var_id,INE,start = (/1,1/), count=(/3, ne_global/))
      CALL GENERIC_NETCDF_ERROR(19, iret)


      iret=nf90_inq_varid(ncid,'ocean_time',varTime_id)
      CALL GENERIC_NETCDF_ERROR(20, iret)
      iret=nf90_inq_varid(ncid,'Hwave',varHs_id)
      CALL GENERIC_NETCDF_ERROR(21, iret)


      ALLOCATE(Aread(np_global))
      ALLOCATE(Bread(np_global))
      ALLOCATE(Cread(np_global))
      FHNDL=10
      open(FHNDL, file=TRIM(BinFile), form='UNFORMATTED')
      idx=0
      DO
        READ(FHNDL,IOSTAT=ISTAT) TIME_4
        IF (ISTAT /= 0) THEN
          EXIT
        ENDIF
        READ(FHNDL) (Aread(IP), Bread(IP), Cread(IP), IP=1,np_global)
        idx=idx+1
        !
        eTime(1)=DBLE(TIME_4)/DAY2SEC
        iret=nf90_put_var(ncid,varTime_id,eTime,start = (/idx/), count=(/1/))
        CALL GENERIC_NETCDF_ERROR(22, iret)
        iret=nf90_put_var(ncid,varHs_id,Cread,start = (/1, idx/), count=(/np_global, 1/))
        CALL GENERIC_NETCDF_ERROR(23, iret)
        !
        avgHs=sum(Cread)/REAL(np_global)
        Print *, 'idx=', idx, 'maxHs=', maxval(Cread), 'avgHs=', avgHs
      END DO
      close(FHNDL)

      iret = NF90_CLOSE(ncid)
      CALL GENERIC_NETCDF_ERROR(24, iret)
      END PROGRAM
!**********************************************************************
!*                                                                    *
!**********************************************************************    
      SUBROUTINE GENERIC_NETCDF_ERROR(idx, iret)
      USE NETCDF
      implicit none
      integer, intent(in) :: iret, idx
      character(len=1000) :: CHRERR
      IF (iret .NE. nf90_noerr) THEN
        CHRERR = nf90_strerror(iret)
        Print *, 'Error at -', idx, '-', TRIM(CHRERR)
        STOP
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************    
