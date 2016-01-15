      MODULE single_netcdf
#ifdef SINGLE_NETCDF_OUTPUT
      USE NETCDF
      USE schism_glbl, ONLY: rkind, errmsg, npa, np, iplg, nvrt, dt
      USE schism_glbl, ONLY: np_global, ne_global, out_rkind
      USE schism_msgp, ONLY : comm, ierr, itype, rtype, myrank, nproc, istatus
      USE schism_glbl, ONLY: iof, ihconsv, isconsv, noutput, indx_out
      IMPLICIT NONE
      integer, allocatable :: ListIPLG(:)
      integer, allocatable :: ListFirst(:)
      integer, allocatable :: ListMNP(:)
      integer, allocatable :: ListNP_RES(:)
      integer, allocatable :: netcdf_his1_rqst(:)
      integer, allocatable :: netcdf_his1_stat(:,:)
      integer, allocatable :: netcdf_his1_type(:)
      integer, allocatable :: netcdf_his2_rqst(:)
      integer, allocatable :: netcdf_his2_stat(:,:)
      integer, allocatable :: netcdf_his2_type(:)
      integer, allocatable :: netcdf_hisN_rqst(:)
      integer, allocatable :: netcdf_hisN_stat(:,:)
      integer, allocatable :: netcdf_hisN_type(:)
      integer NF90_RUNTYPE
      integer NF90_OUTTYPE_HIS
      integer recs_his
      integer istat
      REAL(8) eTimeStart
      !
      ! Timings
      !
      CONTAINS
      SUBROUTINE SCHISM_DATE2JD(year, month, day, hour, min, sec, eJD)
      IMPLICIT NONE
      integer, intent(in) :: year, month, day, hour, min, sec
      real(8), intent(out) :: eJD
      real(8) :: eJDbase, eFracDay
      integer a, y, m
      a = floor((DBLE(14) - DBLE(month))/DBLE(12));
      y = year + 4800 - a;
      m = month + 12*a - 3;
      eJDbase = DBLE(day)                                            &
     & + DBLE(floor((DBLE(153)*DBLE(m) + DBLE(2))/DBLE(5)))          &
     & + DBLE(y)*DBLE(365)                                           &
     & + DBLE(floor(DBLE(y)/DBLE(4)))                                &
     & - DBLE(floor(DBLE(y)/DBLE(100)))                              &
     & + DBLE(floor(DBLE(y)/DBLE(400))) - DBLE(32045)
      eFracDay=(DBLE(sec) +                                          &
     &          DBLE(60)*DBLE(min) +                                 &
     &          DBLE(3600)*(DBLE(hour) - DBLE(12))                   &
     &          )/DBLE(86400)
      eJD=eJDbase + eFracDay
      END SUBROUTINE
      !
      !
      !
      SUBROUTINE SCHISM_ConvertSix2mjd(year, month, day, hour, min, sec, eMJD)
      IMPLICIT NONE
      integer, intent(in) :: year, month, day, hour, min, sec
      real(8), intent(out) :: eMJD
      real(8) :: eJD1, eJD2
      CALL SCHISM_DATE2JD(year, month, day, hour, min, sec, eJD1)
      CALL SCHISM_DATE2JD(1858, 11, 17, 0, 0, 0, eJD2)
      eMJD=eJD1-eJD2
      END SUBROUTINE
      !
      !
      !
      SUBROUTINE SCHISM_ConvertSix2string(year, month, day, hour, min, sec, eTimeStr)
      IMPLICIT NONE
      integer, intent(in) :: year, month, day, hour, min, sec
      character(len=15), intent(out) :: eTimeStr
      WRITE(eTimeStr, 20) year, month, day, hour, min, sec
  20  FORMAT (i4.4, i2.2, i2.2, '.', i2.2, i2.2, i2.2)
      END SUBROUTINE
      !
      !
      !
      SUBROUTINE SCHISM_MONTH_LEN(year, month, lenmonth)
      IMPLICIT NONE
      integer, intent(in) :: year, month
      integer, intent(out) :: lenmonth
      IF ((month .eq. 1).or.(month .eq. 3).or.(month .eq. 5).or.(month .eq. 7).or.(month .eq. 8).or.(month .eq. 10).or.(month .eq. 12)) THEN
        lenmonth=31
      END IF
      IF ((month .eq. 4).or.(month .eq. 6).or.(month .eq. 9).or.(month .eq. 11)) THEN
        lenmonth=30
      END IF
      IF (month .eq. 2) THEN
        IF (MOD(year, 4) .ne. 0) THEN
          lenmonth=28
        ELSE
          IF (MOD(year, 100) .ne. 0) THEN
            lenmonth=29
          ELSE
            IF (MOD(year, 400) .ne. 0) THEN
              lenmonth=28
            ELSE
              lenmonth=29
            END IF
          END IF
        END IF
      END IF
      END SUBROUTINE
      !
      !
      !
      SUBROUTINE SCHISM_JD2DATE(year, month, day, hour, min, sec, eJD)
      ! The following algorithm is from the Calendar FAQ. 
      IMPLICIT NONE
      integer, intent(out) :: year, month, day, hour, min, sec
      real(8), intent(in) :: eJD
      integer ijd, a, b, c, d, e, m
      integer secNear, lenmonth
      real(8) :: fjd, second
      ijd = floor(eJD + 0.5_8)
      a = ijd + 32044
      b = floor((DBLE(4)*DBLE(a) + DBLE(3)) / DBLE(146097))
      c = a - floor((DBLE(b) * DBLE(146097)) / DBLE(4))
      d = floor((DBLE(4)*DBLE(c) + DBLE(3)) / DBLE(1461))
      e = c - floor((DBLE(1461)*DBLE(d)) / DBLE(4));
      m = floor((DBLE(5) * DBLE(e) + DBLE(2)) / DBLE(153))
      day   = e - floor((DBLE(153) * DBLE(m) + DBLE(2)) / DBLE(5)) + 1;
      month = m + 3 - 12 * floor(DBLE(m) / DBLE(10))
      year  = b * 100 + d - 4800 + floor(DBLE(m) / DBLE(10))
      fjd    = eJD - DBLE(ijd) + 0.5_8
      second = DBLE(86400) * fjd
      hour   = floor(second/DBLE(3600))
      second = second - DBLE(3600)*DBLE(hour)
      min    = floor(second/DBLE(60))
      sec    = floor(second - DBLE(60)*min)
      secNear=NINT(second - DBLE(60)*min)
      IF (secNear .eq. 60) THEN
        sec=0
        min=min+1
      END IF
      IF (min .eq. 60) THEN
        min=0
        hour=hour+1
      END IF
      IF (hour .eq. 24) THEN
        hour=0
        day=day+1
      END IF
      CALL SCHISM_MONTH_LEN(year, month, lenmonth)
      IF (day .eq. lenmonth+1) THEN
        day=1
        month=month+1
      END IF
      IF (month .eq. 13) THEN
        month=1
        year=year+1
      END IF
      END SUBROUTINE
      !
      !
      !
      SUBROUTINE SCHISM_MJD2CT(XMJD,STIME)
      IMPLICIT NONE
      CHARACTER(LEN=15), INTENT(OUT) :: STIME
      real(8), INTENT(IN) :: XMJD
      integer year, month, day, hour, min, sec
      real(8) XMJD_1858, eMJD
      CALL SCHISM_DATE2JD(1858, 11, 17, 0, 0, 0, XMJD_1858)
      eMJD = XMJD + XMJD_1858
      CALL SCHISM_JD2DATE(year, month, day, hour, min, sec, eMJD)
      CALL SCHISM_ConvertSix2string(year, month, day, hour, min, sec, STIME)
      END SUBROUTINE
      !
      ! netcdf error
      ! 
      SUBROUTINE GENERIC_NETCDF_ERROR_SCHISM(CallFct, idx, iret)
      implicit none
      integer, intent(in) :: iret, idx
      character(*), intent(in) :: CallFct
      character(len=500) :: CHRERR
      IF (iret .NE. nf90_noerr) THEN
        CHRERR = nf90_strerror(iret)
        WRITE(errmsg,*) 'NETCDF error in routine ', TRIM(CallFct), ' Error Message: ', TRIM(CHRERR), ' Position in the routine :', idx
        CALL parallel_abort(errmsg)
      ENDIF
      END SUBROUTINE
      !
      ! ListFirst and ListIPLG
      !
      SUBROUTINE SCHISM_COLLECT_ALL_IPLG
      IMPLICIT NONE
      integer, allocatable :: rbuf_int(:)
      integer len, iProc, IP, idx, sumMNP
      allocate(ListMNP(nproc), ListNP_RES(nproc), rbuf_int(2), stat=istat)
      IF (istat/=0) CALL parallel_abort('error allocating ListMNP and co')
      IF (myrank == 0) THEN
        ListMNP(1)=npa
        ListNP_RES(1)=np
        DO iProc=2,nproc
          CALL MPI_RECV(rbuf_int,2,itype, iProc-1, 257, comm, istatus, ierr)
          ListMNP(iProc)=rbuf_int(1)
          ListNP_RES(iProc)=rbuf_int(2)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(ListMNP,nproc,itype, iProc-1, 263, comm, ierr)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(ListNP_RES,nproc,itype, iProc-1, 571, comm, ierr)
        END DO
      ELSE
        rbuf_int(1)=npa
        rbuf_int(2)=np
        CALL MPI_SEND(rbuf_int,2,itype, 0, 257, comm, ierr)
        CALL MPI_RECV(ListMNP,nproc,itype, 0, 263, comm, istatus, ierr)
        CALL MPI_RECV(ListNP_RES,nproc,itype, 0, 571, comm, istatus, ierr)
      END IF
      deallocate(rbuf_int)
      sumMNP=sum(ListMNP)
      allocate(ListIPLG(sumMNP), stat=istat)
      IF (istat/=0) CALL parallel_abort('error allocating ListIPLG')
      IF (myrank == 0) THEN
        idx=0
        DO IP=1,npa
          idx=idx+1
          ListIPLG(idx)=iplg(IP)
        END DO
        DO iProc=2,nproc
          len=ListMNP(iProc)
          allocate(rbuf_int(len), stat=istat)
          IF (istat/=0) CALL parallel_abort('error allocating rbuf_int')
          CALL MPI_RECV(rbuf_int,len,itype, iProc-1, 269, comm, istatus, ierr)
          DO IP=1,len
            idx=idx+1
            ListIPLG(idx)=rbuf_int(IP)
          END DO
          deallocate(rbuf_int)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(ListIPLG,sumMNP,itype, iProc-1, 271, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(iplg,npa,itype, 0, 269, comm, ierr)
        CALL MPI_RECV(ListIPLG,sumMNP,itype, 0, 271, comm, istatus, ierr)
      END IF
      ListFirst=0
      DO iProc=2,nproc
         ListFirst(iProc)=ListFirst(iProc-1) + ListMNP(iProc-1)
      END DO
      END SUBROUTINE
      !
      ! 
      !
      SUBROUTINE INIT_PARALLEL_ARRAYS
      IMPLICIT NONE
      integer iProc, IP, IP_glob
      integer eSize, NP_RESloc
      integer, allocatable :: dspl_his1(:), dspl_his2(:), dspl_hisN(:)
      eSize = SIZE(istatus)
      IF (myrank .eq. 0) THEN
        allocate(netcdf_his1_rqst(nproc-1), netcdf_his1_stat(eSize,nproc-1), netcdf_his1_type(nproc-1), stat=istat)
        IF (istat/=0) CALL parallel_abort('Error allocating the netcdf_his1_ arrays')
        allocate(netcdf_his2_rqst(nproc-1), netcdf_his2_stat(eSize,nproc-1), netcdf_his2_type(nproc-1), stat=istat)
        IF (istat/=0) CALL parallel_abort('Error allocating the netcdf_his2_ arrays')
        allocate(netcdf_hisN_rqst(nproc-1), netcdf_hisN_stat(eSize,nproc-1), netcdf_hisN_type(nproc-1), stat=istat)
        IF (istat/=0) CALL parallel_abort('Error allocating the netcdf_hisN_ arrays')
        DO iProc=2,nproc
          NP_RESloc=ListNP_RES(iProc)
          allocate(dspl_his1(NP_RESloc), dspl_his2(NP_RESloc), dspl_hisN(NP_RESloc), stat=istat)
          IF (istat/=0) CALL parallel_abort('Error allocating dspl_his')
          DO IP=1,NP_RESloc
            IP_glob=ListIPLG(IP+ListFirst(iProc))
            dspl_his1(IP)=IP_glob-1
            dspl_his2(IP)=2*(IP_glob-1)
            dspl_hisN(IP)=nvrt*(IP_glob-1)
          END DO
          call mpi_type_create_indexed_block(NP_RESloc,1   ,dspl_his1,rtype,netcdf_his1_type(iProc-1), ierr)
          call mpi_type_commit(netcdf_his1_type(iProc-1), ierr)
          call mpi_type_create_indexed_block(NP_RESloc,2   ,dspl_his2,rtype,netcdf_his2_type(iProc-1), ierr)
          call mpi_type_commit(netcdf_his2_type(iProc-1), ierr)
          call mpi_type_create_indexed_block(NP_RESloc,nvrt,dspl_hisN,rtype,netcdf_hisN_type(iProc-1), ierr)
          call mpi_type_commit(netcdf_hisN_type(iProc-1), ierr)
          deallocate(dspl_his1, dspl_his2, dspl_hisN)
        END DO
      END IF
      END SUBROUTINE
      !
      !
      !
      SUBROUTINE WRITE_1DVAR_SINGLE(ncid, string, VARin)
      IMPLICIT NONE
      integer, intent(in) :: ncid
      character (len = *), parameter :: CallFct = "WRITE_1DVAR_SINGLE"
      character(*), intent(in) :: string
      REAL(rkind), intent(in) :: VARin(npa)
      REAL(rkind), allocatable :: VARout(:)
      integer iProc, IP, IPglob
      integer iret, var_id
      IF (myrank .eq. 0) THEN
        allocate(VARout(np_global), stat=istat)
        DO iProc=2,nproc
          call mpi_irecv(VARout,1,netcdf_his1_type(iProc-1),iProc-1,8024,comm,netcdf_his1_rqst(iProc-1),ierr)
        END DO
        DO IP=1,np
          IPglob = iplg(IP)
          VARout(IPglob) = VARin(IP)
        END DO
        iret=nf90_inq_varid(ncid, TRIM(string), var_id)
        CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 1, iret)
        IF (nproc > 1) THEN
          call mpi_waitall(nproc-1,netcdf_his1_rqst,netcdf_his1_stat,ierr)
        END IF
        IF (NF90_RUNTYPE == NF90_OUTTYPE_HIS) THEN
          iret=nf90_put_var(ncid,var_id,VARout,start = (/1, recs_his/), count = (/ np_global, 1 /))
          CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 2, iret)
        ELSE
          iret=nf90_put_var(ncid,var_id,SNGL(VARout),start = (/1, recs_his/), count = (/ np_global, 1 /))
          CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 3, iret)
        ENDIF
      ELSE
        CALL MPI_SEND(VARin, np, rtype, 0, 8024, comm, ierr)
      END IF
      END SUBROUTINE
      !
      !
      !
      SUBROUTINE DEFINE_1DVAR_SINGLE(ncid, string)
      IMPLICIT NONE
      integer, intent(in) :: ncid
      character(*) :: string
      character (len = *), parameter :: CallFct = "DEFINE_1DVAR_SINGLE"
      integer iret, ntime_dims, mnp_dims, var_id
      iret=nf90_inq_dimid(ncid, 'ocean_time', ntime_dims)
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 1, iret)
      iret=nf90_inq_dimid(ncid, 'np_global', mnp_dims)
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 2, iret)
      iret=nf90_def_var(ncid,TRIM(string),NF90_RUNTYPE,(/ mnp_dims, ntime_dims/), var_id)
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 3, iret)
      END SUBROUTINE
      !
      !
      !
      SUBROUTINE WRITE_NETCDF_TIME_HEADER_SCHISM(ncid, nbTime, ntime_dims)
      IMPLICIT NONE
      integer, intent(in) :: ncid, nbTime
      integer, intent(inout) :: ntime_dims
      character (len = *), parameter :: UNITS = "units"
      character (len = *), parameter :: CallFct="WRITE_NETCDF_TIME_HEADER_SCHISM"
      integer iret, fifteen_dims, var_id
      iret=nf90_inq_dimid(ncid, 'fifteen', fifteen_dims)
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 1, iret)
      IF (nbTime.gt.0) THEN
        iret = nf90_def_dim(ncid, 'ocean_time', nbTime, ntime_dims)
        CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 2, iret)
      ELSE
        iret = nf90_def_dim(ncid, 'ocean_time', NF90_UNLIMITED, ntime_dims)
        CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 3, iret)
      END IF
      iret=nf90_def_var(ncid,'ocean_time',NF90_RUNTYPE,(/ ntime_dims/), var_id)
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 4, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'seconds since 1858-11-17 00:00:00')
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 5, iret)
      iret=nf90_put_att(ncid,var_id,"calendar",'gregorian')
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 6, iret)
      iret=nf90_def_var(ncid,'ocean_time_day',NF90_RUNTYPE,(/ ntime_dims/), var_id)
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 7, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'days since 1858-11-17 00:00:00')
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 8, iret)
      iret=nf90_put_att(ncid,var_id,"calendar",'gregorian')
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 9, iret)
      iret=nf90_def_var(ncid,'ocean_time_str',NF90_CHAR,(/ fifteen_dims, ntime_dims/), var_id)
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 10, iret)
      END SUBROUTINE
      !
      !
      !
      SUBROUTINE WRITE_NETCDF_TIME_SCHISM(ncid, idx, eTimeDay)
      USE NETCDF
      IMPLICIT NONE
      integer, intent(in) :: ncid, idx
      REAL(rkind), intent(IN) :: eTimeDay
      character (len = *), parameter :: CallFct="WRITE_NETCDF_TIME_SCHISM"
      integer oceantimeday_id, oceantimestr_id, oceantime_id
      integer iret, I
      CHARACTER :: eChar
      REAL(rkind) eTimeSec
      CHARACTER(LEN=15) :: eTimeStr
      CALL SCHISM_MJD2CT(eTimeDay,eTimeStr)
      eTimeSec=eTimeDay / 86400.0_rkind
      iret=nf90_inq_varid(ncid, 'ocean_time', oceantime_id)
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 1, iret)
      iret=nf90_put_var(ncid,oceantime_id,eTimeSec,start=(/idx/) )
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 2, iret)
      iret=nf90_inq_varid(ncid, 'ocean_time_day', oceantimeday_id)
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 3, iret)
      iret=nf90_put_var(ncid,oceantimeday_id,eTimeDay,start=(/idx/) )
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 4, iret)
      iret=nf90_inq_varid(ncid, 'ocean_time_str', oceantimestr_id)
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 5, iret)
      DO i=1,15
        eChar=eTimeStr(i:i)
        iret=nf90_put_var(ncid,oceantimestr_id,eChar,start=(/i, idx/) )
        CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 6, iret)
      END DO
      END SUBROUTINE
      !
      !
      !
      SUBROUTINE INIT_SINGLE_OUTPUT
      IMPLICIT NONE
      integer iret, ncid
      character (len = *), parameter :: CallFct="INIT_SINGLE_OUTPUT"
      character(len=*), parameter :: FILE_NAME="outputs/schism_history.nc"
      integer one_dims, two_dims, three_dims, fifteen_dims, mnp_dims, mne_dims
      integer ntime_dims
      integer nbTime
      integer j
      !
      iret = nf90_create(TRIM(FILE_NAME), NF90_CLOBBER, ncid)
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 1, iret)
      iret = nf90_def_dim(ncid, 'one', 1, one_dims)
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 2, iret)
      iret = nf90_def_dim(ncid, 'two', 2, two_dims)
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 3, iret)
      iret = nf90_def_dim(ncid, 'three', 3, three_dims)
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 4, iret)
      iret = nf90_def_dim(ncid, 'fifteen', 15, fifteen_dims)
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 5, iret)
      iret = nf90_def_dim(ncid, 'np_global', np_global, mnp_dims)
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 6, iret)
      iret = nf90_def_dim(ncid, 'ne_global', ne_global, mne_dims)
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 7, iret)
      !
      nbTime=-1
      CALL WRITE_NETCDF_TIME_HEADER_SCHISM(ncid, nbTime, ntime_dims)
      !
      DO j=1,noutput
        if (iof(j)==1) then
          if(j<=13) then
            if(j==1) then
              CALL DEFINE_1DVAR_SINGLE(ncid, "eta2")
            else if(j==2) then
              CALL DEFINE_1DVAR_SINGLE(ncid, "pr")
            else if(j==3.and.ihconsv/=0) then
              CALL DEFINE_1DVAR_SINGLE(ncid, "airt1")
            else if(j==4.and.ihconsv/=0) then
              CALL DEFINE_1DVAR_SINGLE(ncid, "shum1")
            else if(j==5.and.ihconsv/=0) then
              CALL DEFINE_1DVAR_SINGLE(ncid, "srad")
            else if(j==6.and.ihconsv/=0) then
              CALL DEFINE_1DVAR_SINGLE(ncid, "fluxsu")
            else if(j==7.and.ihconsv/=0) then
              CALL DEFINE_1DVAR_SINGLE(ncid, "fluxlu")
            else if(j==8.and.ihconsv/=0) then
              CALL DEFINE_1DVAR_SINGLE(ncid, "hradu")
            else if(j==9.and.ihconsv/=0) then
              CALL DEFINE_1DVAR_SINGLE(ncid, "hradd")
            else if(j==10.and.ihconsv/=0) then
              CALL DEFINE_1DVAR_SINGLE(ncid, "sflux")
            else if(j==11.and.isconsv/=0) then
              CALL DEFINE_1DVAR_SINGLE(ncid, "fluxevp")
            else if(j==12.and.isconsv/=0) then
              CALL DEFINE_1DVAR_SINGLE(ncid, "fluxprc")
            else if(j==13) then
              CALL DEFINE_1DVAR_SINGLE(ncid, "Cdp")
            endif
          else if(j<=16) then
            if(j==14) then
              CALL DEFINE_1DVAR_SINGLE(ncid, "windx")
              CALL DEFINE_1DVAR_SINGLE(ncid, "windy")
            else if(j==15) then !in ll frame if ics=2
              CALL DEFINE_1DVAR_SINGLE(ncid, "taux")
              CALL DEFINE_1DVAR_SINGLE(ncid, "tauy")
            else !j=16
              CALL DEFINE_1DVAR_SINGLE(ncid, "davx")
              CALL DEFINE_1DVAR_SINGLE(ncid, "davy")
            endif
          else if(j<27) then
            call parallel_abort('No code for netcdf vertical arrays yet 1')
          else if(j==27) then
            call parallel_abort('No code for netcdf vertical arrays yet 2')
          else !optional modules; MUST BE IN THE SAME ORDER AS BEFORE

#ifdef USE_GEN
            if(j>=indx_out(1,1).and.j<=indx_out(1,2)) then
              call parallel_abort('No code for netcdf vertical arrays yet 3')
            endif
#endif

#ifdef USE_AGE
            if(j>=indx_out(2,1).and.j<=indx_out(2,2)) then
              call parallel_abort('No code for netcdf vertical arrays yet 4')
            endif
#endif

#ifdef USE_SED
            if(j>=indx_out(3,1).and.j<=indx_out(3,2)) then
              call parallel_abort('No code for netcdf sediment yet')
            endif !scope of SED model
#endif /*USE_SED*/

#ifdef USE_ECO
            if(j>=indx_out(4,1).and.j<=indx_out(4,2)) then
              call parallel_abort('No code for netcdf ECO yet')
            endif
#endif

#ifdef USE_ICM
            if(j>=indx_out(5,1).and.j<=indx_out(5,2)) then
              call parallel_abort('No code for netcdf ICM yet')
            endif
#endif

#ifdef USE_COSINE
            if(j>=indx_out(8,1).and.j<=indx_out(8,2)) then
              call parallel_abort('No code for netcdf COSINE yet')
            endif
#endif

#ifdef USE_FIB
            if(j>=indx_out(9,1).and.j<=indx_out(9,2)) then
              call parallel_abort('No code for netcdf FIB yet')
            endif
#endif

#ifdef USE_TIMOR
            if(j>=indx_out(10,1).and.j<=indx_out(10,2)) then
              call parallel_abort('No code for netcdf TIMOR yet')
            endif
#endif

#ifdef USE_SED2D
            if(j>=indx_out(6,1).and.j<=indx_out(6,2)) then
              call parallel_abort('No code for netcdf SED2D yet')
            endif !scope of SED2D model
#endif /*USE_SED2D*/

#ifdef USE_WWM
            if((j>=indx_out(7,1)).and.(j<=indx_out(7,2))) then
              call parallel_abort('No code for netcdf WWM yet')
            endif !scope of WWM; j<=indx_out(3,2)
#endif /*USE_WWM*/
          end if
        end if
      END DO
      END SUBROUTINE
      !
      !
      !
      SUBROUTINE WRITE_SINGLE_OUTPUT_DATA(it)
      USE schism_glbl, only : eta2, pr, airt1, shum1, srad, fluxsu, fluxlu, hradu, hradd
      USE schism_glbl, only : sflux, fluxevp, fluxprc, Cdp, windx, windy, tau, dav
      IMPLICIT NONE
      integer, intent(in) :: it
      integer iret, ncid
      character (len = *), parameter :: CallFct="WRITE_SINGLE_OUTPUT_DATA"
      character(len=*), parameter :: FILE_NAME='outputs/schism_history.nc'
      real(rkind) :: eTimeDay
      integer j

      !
      iret=nf90_open(TRIM(FILE_NAME), nf90_write, ncid)
      CALL GENERIC_NETCDF_ERROR_SCHISM(CallFct, 1, iret)
      !
      eTimeDay = eTimeStart + it * (dt/86400.)
      CALL WRITE_NETCDF_TIME_SCHISM(ncid, recs_his, eTimeDay)
      !
      DO j=1,noutput
        if (iof(j)==1) then
          if(j<=13) then
            if(j==1) then
              CALL WRITE_1DVAR_SINGLE(ncid, "eta2", eta2)
            else if(j==2) then
              CALL WRITE_1DVAR_SINGLE(ncid, "pr", pr)
            else if(j==3.and.ihconsv/=0) then
              CALL WRITE_1DVAR_SINGLE(ncid, "airt1", airt1)
            else if(j==4.and.ihconsv/=0) then
              CALL WRITE_1DVAR_SINGLE(ncid, "shum1", shum1)
            else if(j==5.and.ihconsv/=0) then
              CALL WRITE_1DVAR_SINGLE(ncid, "srad", srad)
            else if(j==6.and.ihconsv/=0) then
              CALL WRITE_1DVAR_SINGLE(ncid, "fluxsu", fluxsu)
            else if(j==7.and.ihconsv/=0) then
              CALL WRITE_1DVAR_SINGLE(ncid, "fluxlu", fluxlu)
            else if(j==8.and.ihconsv/=0) then
              CALL WRITE_1DVAR_SINGLE(ncid, "hradu", hradu)
            else if(j==9.and.ihconsv/=0) then
              CALL WRITE_1DVAR_SINGLE(ncid, "hradd", hradd)
            else if(j==10.and.ihconsv/=0) then
              CALL WRITE_1DVAR_SINGLE(ncid, "sflux", sflux)
            else if(j==11.and.isconsv/=0) then
              CALL WRITE_1DVAR_SINGLE(ncid, "fluxevp", fluxevp)
            else if(j==12.and.isconsv/=0) then
              CALL WRITE_1DVAR_SINGLE(ncid, "fluxprc", fluxprc)
            else if(j==13) then
              CALL WRITE_1DVAR_SINGLE(ncid, "Cdp", Cdp)
            endif
          else if(j<=16) then
            if(j==14) then
              CALL WRITE_1DVAR_SINGLE(ncid, "windx", windx)
              CALL WRITE_1DVAR_SINGLE(ncid, "windy", windy)
            else if(j==15) then !in ll frame if ics=2
              CALL WRITE_1DVAR_SINGLE(ncid, "taux", tau(1,:))
              CALL WRITE_1DVAR_SINGLE(ncid, "tauy", tau(2,:))
            else !j=16
              CALL WRITE_1DVAR_SINGLE(ncid, "davx", dav(1,:))
              CALL WRITE_1DVAR_SINGLE(ncid, "davy", dav(2,:))
            endif
          else if(j<27) then
            call parallel_abort('No code for netcdf vertical arrays yet 1')
          else if(j==27) then
            call parallel_abort('No code for netcdf vertical arrays yet 2')
          else !optional modules; MUST BE IN THE SAME ORDER AS BEFORE

#ifdef USE_GEN
            if(j>=indx_out(1,1).and.j<=indx_out(1,2)) then
              call parallel_abort('No code for netcdf vertical arrays yet 3')
            endif
#endif

#ifdef USE_AGE
            if(j>=indx_out(2,1).and.j<=indx_out(2,2)) then
              call parallel_abort('No code for netcdf vertical arrays yet 4')
            endif
#endif

#ifdef USE_SED
            if(j>=indx_out(3,1).and.j<=indx_out(3,2)) then
              call parallel_abort('No code for netcdf sediment yet')
            endif !scope of SED model
#endif /*USE_SED*/

#ifdef USE_ECO
            if(j>=indx_out(4,1).and.j<=indx_out(4,2)) then
              call parallel_abort('No code for netcdf ECO yet')
            endif
#endif

#ifdef USE_ICM
            if(j>=indx_out(5,1).and.j<=indx_out(5,2)) then
              call parallel_abort('No code for netcdf ICM yet')
            endif
#endif

#ifdef USE_COSINE
            if(j>=indx_out(8,1).and.j<=indx_out(8,2)) then
              call parallel_abort('No code for netcdf COSINE yet')
            endif
#endif

#ifdef USE_FIB
            if(j>=indx_out(9,1).and.j<=indx_out(9,2)) then
              call parallel_abort('No code for netcdf FIB yet')
            endif
#endif

#ifdef USE_TIMOR
            if(j>=indx_out(10,1).and.j<=indx_out(10,2)) then
              call parallel_abort('No code for netcdf TIMOR yet')
            endif
#endif

#ifdef USE_SED2D
            if(j>=indx_out(6,1).and.j<=indx_out(6,2)) then
              call parallel_abort('No code for netcdf SED2D yet')
            endif !scope of SED2D model
#endif /*USE_SED2D*/

#ifdef USE_WWM
            if((j>=indx_out(7,1)).and.(j<=indx_out(7,2))) then
              call parallel_abort('No code for netcdf WWM yet')
            endif !scope of WWM; j<=indx_out(3,2)
#endif /*USE_WWM*/
          end if
        end if
      END DO
      END SUBROUTINE
      !
      ! 
      !
      SUBROUTINE GENERAL_INIT_NETCDF_OUTPUT_SCHISM
      IMPLICIT NONE
      REAL(8) eMJD
      integer :: sim_year = 0
      integer :: sim_month = 0
      integer :: sim_day = 0
      integer :: sim_hour = 0
      integer :: sim_minute = 0
      integer :: sim_second = 0
      IF (out_rkind == 4) THEN
        NF90_OUTTYPE_HIS = NF90_REAL
      ELSE
        NF90_OUTTYPE_HIS = NF90_DOUBLE
      END IF
      IF (out_rkind == 4) THEN
        NF90_RUNTYPE = NF90_REAL
      ELSE
        NF90_RUNTYPE = NF90_DOUBLE
      END IF
      recs_his = 0
      CALL SCHISM_COLLECT_ALL_IPLG
      CALL INIT_PARALLEL_ARRAYS
      CALL INIT_SINGLE_OUTPUT
      CALL SCHISM_ConvertSix2mjd(sim_year, sim_month, sim_day, sim_hour, sim_minute, sim_second, eMJD)
      IF (rkind .eq. 4) THEN
        eTimeStart = SNGL(eMJD)
      ELSE
        eTimeStart = eMJD
      END IF
      END SUBROUTINE
      !
      ! 
      !
      SUBROUTINE NETCDF_SINGLE_OUTPUT(it)
      IMPLICIT NONE
      integer, intent(in) :: it
      logical :: IsInitialized = .FALSE.
      IF (IsInitialized .eqv. .FALSE.) THEN
        CALL GENERAL_INIT_NETCDF_OUTPUT_SCHISM
        IsInitialized=.TRUE.
      END IF
      recs_his = recs_his + 1
      CALL WRITE_SINGLE_OUTPUT_DATA(it)
      END SUBROUTINE
      !
      ! 
      !
#endif
      END MODULE
