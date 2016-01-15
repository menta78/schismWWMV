!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!**********************************************************************
!*  Author: Ivica Janeković [ivica.jan@gmail.com]                                                                   *
!**********************************************************************
! This is set of subroutines for handling wind and mslp from netcdf WRF
! to get Sschism forcing fields.
! list of SUBROUTINES:
! INIT_NETCDF_CF 	loads wind_time, compute interp coefs.
! READ_INTERP_NETCDF_CF	reads fields and interpolate them on the FEM
! READ_NETCDF_DIRECT    fileds already interpolated into SCHICM grid
! FIX_COORDS      	tranform coords into Sschism radians
!*********************************************************************
      SUBROUTINE INIT_NETCDF_CF
      USE NETCDF
      USE schism_glbl, only : npa, rkind, pi
      USE schism_glbl, only : cf_c11, cf_c21, cf_c22, cf_c12
      USE schism_glbl, only : cf_a, cf_b, cf_c, cf_d, cf_J
      USE schism_glbl, only : cf_add_offset_uwind, cf_scale_factor_uwind
      USE schism_glbl, only : cf_add_offset_vwind, cf_scale_factor_vwind
      USE schism_glbl, only : cf_add_offset_pr, cf_scale_factor_pr
      USE schism_glbl, only : wind_time_sec, nwtimes 
      USE schism_glbl, only : xlon, ylat
      USE schism_glbl, only : NDX_WIND_FD, NDY_WIND_FD
      USE schism_msgp, only : myrank

! hardcoded WIND MSLP file name to UVP.nc !

      IMPLICIT NONE
      integer                  :: I, J, ISTAT, fid, varid, dimids(2), dimidsB(3)
      integer		       :: i11, i12, i21, j11, j12, j21, WindTimeToSec
      real(rkind), allocatable :: dist(:,:)
      integer     closest(2)
      real(rkind) d_lon, d_lat
      REAL(rkind), allocatable :: CF_LON(:,:), CF_LAT(:,:), wind_time_raw(:)
      character (len = *), parameter :: CallFct="INIT_NETCDF_CF"
      character(len=100)       :: CHRERR, WindTimeUnits

      ! Open NC file
      ISTAT = nf90_open('UVP.nc', nf90_nowrite, fid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

      ! Reading wind_time
      ISTAT = nf90_inq_varid(fid, "wind_time", varid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)
      ISTAT = nf90_inquire_variable(fid, varid, dimids=dimids)
      CALL GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)
      ISTAT = nf90_inquire_dimension(fid, dimids(1), len=nwtimes)
      CALL GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)
      allocate(wind_time_raw(nwtimes), stat=istat)
      ISTAT = nf90_get_var(fid, varid, wind_time_raw)
      ISTAT = nf90_get_att(fid, varid, "units", WindTimeUnits)
      CALL CF_EXTRACT_TIME_NEW(WindTimeUnits, WindTimeToSec)
      allocate(wind_time_sec(nwtimes), stat=istat)
      wind_time_sec(:)=(wind_time_raw(:)-wind_time_raw(1))*WindTimeToSec
      deallocate(wind_time_raw)

      ! Reading Uwind attributes
      ISTAT = nf90_inq_varid(fid, "Uwind", varid)
      IF (ISTAT /= 0) THEN
      CHRERR = nf90_strerror(ISTAT)
      if(myrank==0)  WRITE(16,*) 'There is no Uwind variable=', TRIM(CHRERR)
      ENDIF
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)
      ! scale_factor
      ISTAT = nf90_get_att(fid, varid, "scale_factor", cf_scale_factor_uwind)
      IF (ISTAT /= 0) THEN
        CHRERR = nf90_strerror(ISTAT)
	cf_scale_factor_uwind=1
      if(myrank==0)  WRITE(16,*) 'CHRERR=', TRIM(CHRERR)
      ENDIF
      if(myrank==0) WRITE(16,*) 'cf_scale_factor_uwind=', cf_scale_factor_uwind
      ! add_offset
      ISTAT = nf90_get_att(fid, varid, "add_offset", cf_add_offset_uwind)
      IF (ISTAT /= 0) THEN
	cf_add_offset_uwind=0
        CHRERR = nf90_strerror(ISTAT)
      if(myrank==0)   WRITE(16,*) 'CHRERR=', TRIM(CHRERR)
      ENDIF
      if(myrank==0) WRITE(16,*) 'cf_add_offset_uwind=', cf_add_offset_uwind

      ! Reading Vwind attributes
      ISTAT = nf90_inq_varid(fid, "Vwind", varid)
      IF (ISTAT /= 0) THEN
        CHRERR = nf90_strerror(ISTAT)
      if(myrank==0)  WRITE(16,*) 'There is no Vwind variable=', TRIM(CHRERR)
      ENDIF
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)
      ISTAT = nf90_get_att(fid, varid, "scale_factor", cf_scale_factor_vwind)
      IF (ISTAT /= 0) THEN
        CHRERR = nf90_strerror(ISTAT)
        if(myrank==0) WRITE(16,*) 'CHRERR=', TRIM(CHRERR)
        cf_scale_factor_vwind=1
      ENDIF
      if(myrank==0) WRITE(16,*) 'cf_scale_factor_vwind=', cf_scale_factor_vwind

      ISTAT = nf90_get_att(fid, varid, "add_offset", cf_add_offset_vwind)
      IF (ISTAT /= 0) THEN
        CHRERR = nf90_strerror(ISTAT)
        if(myrank==0) WRITE(16,*) 'CHRERR=', TRIM(CHRERR)
        cf_add_offset_vwind=0
      ENDIF
      if(myrank==0) WRITE(16,*) 'cf_add_offset_vwind=', cf_add_offset_vwind

      ! Reading Pair attributes
      ISTAT = nf90_inq_varid(fid, "Pair", varid)
      IF (ISTAT /= 0) THEN
      CHRERR = nf90_strerror(ISTAT)
      if(myrank==0)  WRITE(16,*) 'There is no Pair variable !!!', TRIM(CHRERR)
      ENDIF
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)
      ISTAT = nf90_get_att(fid, varid, "scale_factor", cf_scale_factor_pr)
      IF (ISTAT /= 0) THEN
        CHRERR = nf90_strerror(ISTAT)
        if(myrank==0) WRITE(16,*) 'CHRERR=', TRIM(CHRERR)
        cf_scale_factor_pr=1
      ENDIF
      if(myrank==0) WRITE(16,*) 'cf_scale_factor_pr=', cf_scale_factor_pr
      ISTAT = nf90_get_att(fid, varid, "add_offset", cf_add_offset_pr)
      IF (ISTAT /= 0) THEN
        CHRERR = nf90_strerror(ISTAT)
        if(myrank==0) WRITE(16,*) 'CHRERR=', TRIM(CHRERR)
        cf_add_offset_pr=0
      ENDIF
      if(myrank==0) WRITE(16,*) 'cf_add_offset_pr=', cf_add_offset_pr
      
     ! Reading lontitude/latitude array dimensions only
      ISTAT = nf90_inq_varid(fid, "lon2d", varid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)
      ISTAT = nf90_inquire_variable(fid, varid, dimids=dimids)
      CALL GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)
      ISTAT = nf90_inquire_dimension(fid, dimids(1), len=NDX_WIND_FD)
      CALL GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)
      ISTAT = nf90_inquire_dimension(fid, dimids(2), len=NDY_WIND_FD)
      CALL GENERIC_NETCDF_ERROR(CallFct, 9, ISTAT)
!    ALLOCATE METEO COORD ARRAYS ONLY, NO NEED FOR FIELDS LIKE UWIND,VWIND,MSLP
      allocate(CF_LON(NDX_WIND_FD, NDY_WIND_FD), CF_LAT(NDX_WIND_FD, NDY_WIND_FD), stat=istat)

     ! Reading lontitude/latitude arrays
      ISTAT = nf90_inq_varid(fid, "lat2d", varid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 10, ISTAT)
      ISTAT = nf90_get_var(fid, varid, CF_LAT)
      CALL GENERIC_NETCDF_ERROR(CallFct, 11, ISTAT)
      ISTAT = nf90_inq_varid(fid, "lon2d", varid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 12, ISTAT)
      ISTAT = nf90_get_var(fid, varid, CF_LON)
      CALL GENERIC_NETCDF_ERROR(CallFct, 13, ISTAT)
      ISTAT = nf90_close(fid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 19, ISTAT)

      CF_LON(:,:)=CF_LON(:,:)*pi/180.0_rkind
      CF_LAT(:,:)=CF_LAT(:,:)*pi/180.0_rkind


      if(myrank==0) then
      WRITE(16,*) 'nwtimes=', nwtimes
      WRITE(16,*) 'NDX_WIND_FD=', NDX_WIND_FD
      WRITE(16,*) 'NYX_WIND_FD=', NDY_WIND_FD
      WRITE(16,*) 'Min/Max CF_LON',minval(CF_LON),maxval(CF_LON)
      WRITE(16,*) 'Min/Max CF_LAT',minval(CF_LAT),maxval(CF_LAT)
      WRITE(16,*) 'Min/Max xlon',minval(xlon),maxval(xlon)
      WRITE(16,*) 'Min/Max ylat',minval(ylat),maxval(ylat)
      WRITE(16,*) 'Done with UVP.nc init phase, calculating interp coefs'
      FLUSH(16)
      endif

 
! compute nodes and coefs for bilinear interpolation for whole grid
        ALLOCATE(cf_c11(npa,2), cf_c12(npa,2), cf_c21(npa,2), cf_c22(npa,2), stat=istat)
        IF (istat/=0) WRITE(16,*) 'Problem with allocate cf11,12,21,22'
        ALLOCATE(cf_a(npa), cf_b(npa), cf_c(npa), cf_d(npa), cf_J(npa), stat=istat)
        IF (istat/=0) WRITE(16,*) 'Problem with allocate cf_a,b,c,f,J'
        ALLOCATE(dist(NDX_WIND_FD, NDY_WIND_FD), stat=istat)
        IF (istat/=0) WRITE(16,*) 'Problem with allocate dist'

        DO I = 1, npa
          dist(:,:) = ABS( CMPLX(xlon(I)-CF_LON(:,:), ylat(I)-CF_LAT(:,:)) )
          closest(1:2) = MINLOC(dist)
          d_lon = xlon(I)-CF_LON(closest(1),closest(2)) 
          d_lat = ylat(I)-CF_LAT(closest(1),closest(2))
          IF ((d_lon.ge.0).and.(d_lat.ge.0)) THEN ! point is in the I kvadrant
            cf_c11(I,:) = closest(:)
            cf_c21(I,1) = closest(1) + 1
            cf_c22(I,1) = closest(1) + 1
            cf_c12(I,1) = closest(1)
            cf_c21(I,2) = closest(2)
            cf_c22(I,2) = closest(2) + 1
            cf_c12(I,2) = closest(2) + 1
          END IF
          IF ((d_lon.ge.0).and.(d_lat.le.0)) THEN ! point is in the IV kvadrant
            cf_c11(I,1) = closest(1)
            cf_c21(I,1) = closest(1) + 1
            cf_c22(I,1) = closest(1) + 1
            cf_c12(I,:) = closest(:)
            cf_c11(I,2) = closest(2) - 1
            cf_c21(I,2) = closest(2) - 1
            cf_c22(I,2) = closest(2) 
          END IF
          IF ((d_lon.le.0).and.(d_lat.ge.0)) THEN ! point is in the II kvadrant
            cf_c11(I,1) = closest(1) - 1 
            cf_c21(I,:) = closest(:)
            cf_c22(I,1) = closest(1)
            cf_c12(I,1) = closest(1) - 1
            cf_c11(I,2) = closest(2)
            cf_c22(I,2) = closest(2) + 1
            cf_c12(I,2) = closest(2) + 1 
          END IF
          IF ((d_lon.le.0).and.(d_lat.le.0)) THEN ! point is in the III kvadrant
            cf_c11(I,1) = closest(1) - 1
            cf_c21(I,1) = closest(1)
            cf_c22(I,:) = closest(:)
            cf_c12(I,1) = closest(1) - 1
            cf_c11(I,2) = closest(2) - 1
            cf_c21(I,2) = closest(2) - 1
            cf_c12(I,2) = closest(2) 
          END IF
          ! J =1/((x2-x1)*(y2-y1))
          i11=cf_c11(I,1)
          j11=cf_c11(I,2)
          i12=cf_c12(I,1)
          j12=cf_c12(I,2)
          i21=cf_c21(I,1)
          j21=cf_c21(I,2)
          cf_J(I)=1.0/( (CF_LON(i21,j21)-CF_LON(i11,j11))*(CF_LAT(i12,j12)-CF_LAT(i11,j11)) )
          cf_a(I) = CF_LON(i21,j21) - xlon(I) ! x2-x
          cf_b(I) = xlon(I) - CF_LON(i11,j11) ! x-x1
          cf_c(I) = CF_LAT(i12,j12) - ylat(I) ! y2-y
          cf_d(I) = ylat(I) - CF_LAT(i11,j11) ! y-y1
        END DO
      DEALLOCATE(dist)
      DEALLOCATE(CF_LON, CF_LAT)
      if(myrank==0) then
      WRITE(16,*) 'Done with UVP.nc interp coefs, all done in INIT_NETCDF_CF'
      FLUSH(16)
      endif
      END SUBROUTINE INIT_NETCDF_CF
! ***************************************************************************************************
      SUBROUTINE INIT_NETCDF_DIRECT
      USE NETCDF
      USE schism_glbl, only : cf_add_offset_uwind, cf_scale_factor_uwind
      USE schism_glbl, only : cf_add_offset_vwind, cf_scale_factor_vwind
      USE schism_glbl, only : cf_add_offset_pr, cf_scale_factor_pr
      USE schism_glbl, only : wind_time_sec, nwtimes, rkind
      USE schism_msgp, only : myrank

! hardcoded WIND MSLP file name to UVP_direct.nc !

      IMPLICIT NONE
      integer                  :: ISTAT, fid, varid, dimids(2), WindTimeToSec
      REAL(rkind),ALLOCATABLE  :: wind_time_raw(:)
      character (len = *), parameter :: CallFct="INIT_NETCDF_DIRECT"
      character(len=100)       :: CHRERR, WindTimeUnits

      ! Open NC file
      ISTAT = nf90_open('UVP_direct.nc', nf90_nowrite, fid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

     ! Reading wind_time
      ISTAT = nf90_inq_varid(fid, "wind_time", varid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)
      ISTAT = nf90_inquire_variable(fid, varid, dimids=dimids)
      CALL GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)
      ISTAT = nf90_inquire_dimension(fid, dimids(1), len=nwtimes)
      CALL GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)
      allocate(wind_time_raw(nwtimes), stat=istat)
      ISTAT = nf90_get_var(fid, varid, wind_time_raw)
      ISTAT = nf90_get_att(fid, varid, "units", WindTimeUnits)
      CALL CF_EXTRACT_TIME_NEW(WindTimeUnits, WindTimeToSec)
      allocate(wind_time_sec(nwtimes), stat=istat)
      wind_time_sec(:)=(wind_time_raw(:)-wind_time_raw(1))*WindTimeToSec
      deallocate(wind_time_raw)
      
      ! Reading Uwind attributes
      ISTAT = nf90_inq_varid(fid, "Uwind", varid)
      IF (ISTAT /= 0) THEN
      CHRERR = nf90_strerror(ISTAT)
      if(myrank==0)  WRITE(16,*) 'There is no Uwind variable=', TRIM(CHRERR)
      ENDIF
      
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)
      ! scale_factor
      ISTAT = nf90_get_att(fid, varid, "scale_factor", cf_scale_factor_uwind)
      IF (ISTAT /= 0) THEN
        CHRERR = nf90_strerror(ISTAT)
	cf_scale_factor_uwind=1
      if(myrank==0)  WRITE(16,*) 'CHRERR=', TRIM(CHRERR)
      ENDIF
      if(myrank==0) WRITE(16,*) 'cf_scale_factor_uwind=', cf_scale_factor_uwind
      
      ! add_offset
      ISTAT = nf90_get_att(fid, varid, "add_offset", cf_add_offset_uwind)
      IF (ISTAT /= 0) THEN
	cf_add_offset_uwind=0
        CHRERR = nf90_strerror(ISTAT)
      if(myrank==0)   WRITE(16,*) 'CHRERR=', TRIM(CHRERR)
      ENDIF
      if(myrank==0) WRITE(16,*) 'cf_add_offset_uwind=', cf_add_offset_uwind

      ! Reading Vwind attributes
      ISTAT = nf90_inq_varid(fid, "Vwind", varid)
      IF (ISTAT /= 0) THEN
        CHRERR = nf90_strerror(ISTAT)
      if(myrank==0)  WRITE(16,*) 'There is no Vwind variable=', TRIM(CHRERR)
      ENDIF
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)
      ISTAT = nf90_get_att(fid, varid, "scale_factor", cf_scale_factor_vwind)
      IF (ISTAT /= 0) THEN
        CHRERR = nf90_strerror(ISTAT)
        if(myrank==0) WRITE(16,*) 'CHRERR=', TRIM(CHRERR)
        cf_scale_factor_vwind=1
      ENDIF
      if(myrank==0) WRITE(16,*) 'cf_scale_factor_vwind=', cf_scale_factor_vwind

      ISTAT = nf90_get_att(fid, varid, "add_offset", cf_add_offset_vwind)
      IF (ISTAT /= 0) THEN
        CHRERR = nf90_strerror(ISTAT)
        if(myrank==0) WRITE(16,*) 'CHRERR=', TRIM(CHRERR)
        cf_add_offset_vwind=0
      ENDIF
      if(myrank==0) WRITE(16,*) 'cf_add_offset_vwind=', cf_add_offset_vwind

      ! Reading Pair attributes
      ISTAT = nf90_inq_varid(fid, "Pair", varid)
      IF (ISTAT /= 0) THEN
      CHRERR = nf90_strerror(ISTAT)
      if(myrank==0)  WRITE(16,*) 'There is no Pair variable !!!', TRIM(CHRERR)
      ENDIF
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)
      ISTAT = nf90_get_att(fid, varid, "scale_factor", cf_scale_factor_pr)
      IF (ISTAT /= 0) THEN
        CHRERR = nf90_strerror(ISTAT)
        if(myrank==0) WRITE(16,*) 'CHRERR=', TRIM(CHRERR)
        cf_scale_factor_pr=1
      ENDIF
      if(myrank==0) WRITE(16,*) 'cf_scale_factor_pr=', cf_scale_factor_pr

      ISTAT = nf90_get_att(fid, varid, "add_offset", cf_add_offset_pr)
      IF (ISTAT /= 0) THEN
        CHRERR = nf90_strerror(ISTAT)
        if(myrank==0) WRITE(16,*) 'CHRERR=', TRIM(CHRERR)
        cf_add_offset_pr=0
      ENDIF
      if(myrank==0) WRITE(16,*) 'cf_add_offset_pr=', cf_add_offset_pr
 
      if(myrank==0) then
      WRITE(16,*) '... Done with UVP_direct.nc init phase'
      FLUSH(16)
      endif
      END SUBROUTINE INIT_NETCDF_DIRECT
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_INTERP_NETCDF_CF(RECORD_IN,RECORD_OUT)
      USE NETCDF
      USE schism_glbl, only : npa,cf_c11, cf_c21, rkind
      USE schism_glbl, only : cf_c22, cf_c12, cf_a, cf_b, cf_c, cf_d, cf_J
      USE schism_glbl, only : cf_add_offset_uwind, cf_scale_factor_uwind
      USE schism_glbl, only : cf_add_offset_vwind, cf_scale_factor_vwind
      USE schism_glbl, only : cf_add_offset_pr, cf_scale_factor_pr
      USE schism_glbl, only : NDX_WIND_FD, NDY_WIND_FD
      USE schism_msgp, only : myrank

      IMPLICIT NONE
      INTEGER, INTENT(in)                      :: RECORD_IN	! time in [s] from the start of simulation
      REAL(rkind),dimension(npa,3),INTENT(out) :: RECORD_OUT
      REAL(rkind), ALLOCATABLE :: wind_time(:)
      INTEGER                             :: FID, ID, ISTAT, I, dimids(2), ntimes
      REAL(rkind)			  :: UWIND_FD(NDX_WIND_FD, NDY_WIND_FD) 
      REAL(rkind)                         :: VWIND_FD(NDX_WIND_FD, NDY_WIND_FD) 
      REAL(rkind)                         :: MSLP_FD(NDX_WIND_FD, NDY_WIND_FD)
      character (len = *), parameter :: CallFct="READ_INTERP_NETCDF_CF"
      character(len=100)       :: CHRERR

! OPEN FORCING FILE
      ISTAT = NF90_OPEN('UVP.nc', NF90_NOWRITE, FID)
      IF (ISTAT /= 0) THEN
      CHRERR = nf90_strerror(ISTAT)
      if(myrank==0)  WRITE(16,*) 'There is no UVP.nc file!! ', TRIM(CHRERR)
      ENDIF
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

      ISTAT = NF90_inq_varid(FID, 'Uwind', ID)
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)
      IF (ISTAT /= 0) THEN
      CHRERR = nf90_strerror(ISTAT)
      if(myrank==0)  WRITE(16,*) 'There is no Uwind variable ', TRIM(CHRERR)
      UWIND_FD(:,:)=0
      ELSE
      ISTAT = NF90_GET_VAR(FID, ID, UWIND_FD, start = (/ 1, 1, RECORD_IN /), count = (/ NDX_WIND_FD, NDY_WIND_FD, 1 /))
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)
      UWIND_FD(:,:)=cf_add_offset_uwind + cf_scale_factor_uwind*UWIND_FD(:,:)
      ENDIF

      ISTAT = NF90_inq_varid(FID, 'Vwind', ID)
      CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)
      IF (ISTAT /= 0) THEN
      CHRERR = nf90_strerror(ISTAT)
      if(myrank==0)  WRITE(16,*) 'There is no Vwind variable ', TRIM(CHRERR)
      VWIND_FD(:,:)=0
      ELSE
      ISTAT = NF90_GET_VAR(FID, ID, VWIND_FD, start = (/ 1, 1, RECORD_IN /), count = (/ NDX_WIND_FD, NDY_WIND_FD, 1 /))
      CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)
      VWIND_FD(:,:)=cf_add_offset_vwind + cf_scale_factor_vwind*VWIND_FD(:,:)
      ENDIF

      ISTAT = NF90_inq_varid(FID, 'Pair', ID)
      CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)
      IF (ISTAT /= 0) THEN
      CHRERR = nf90_strerror(ISTAT)
      if(myrank==0)  WRITE(16,*) 'There is no Pair variable ', TRIM(CHRERR)
      MSLP_FD(:,:)=101325
      ELSE
      ISTAT = NF90_GET_VAR(FID, ID, MSLP_FD, start = (/ 1, 1, RECORD_IN /), count = (/ NDX_WIND_FD, NDY_WIND_FD, 1 /))
      CALL GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)
      ! MSLP in the WRF file is in mbar but have to use Pa -> multiply with 100.0
      MSLP_FD(:,:)=100*(cf_add_offset_pr+cf_scale_factor_pr*MSLP_FD(:,:))
      ENDIF

      ISTAT = NF90_CLOSE(FID)
      CALL GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)

        DO I = 1, npa
          RECORD_OUT(I,1) =  cf_J(I)*(                                &
     &      UWIND_FD(cf_c11(I,1),cf_c11(I,2))*cf_a(I)*cf_c(I)+        &
     &      UWIND_FD(cf_c21(I,1),cf_c21(I,2))*cf_b(I)*cf_c(I)+        &
     &      UWIND_FD(cf_c12(I,1),cf_c12(I,2))*cf_a(I)*cf_d(I)+        &
     &      UWIND_FD(cf_c22(I,1),cf_c22(I,2))*cf_b(I)*cf_d(I) )
          RECORD_OUT(I,2) = cf_J(I)*(                                 &
     &      VWIND_FD(cf_c11(I,1),cf_c11(I,2))*cf_a(I)*cf_c(I)+        &
     &      VWIND_FD(cf_c21(I,1),cf_c21(I,2))*cf_b(I)*cf_c(I)+        &
     &      VWIND_FD(cf_c12(I,1),cf_c12(I,2))*cf_a(I)*cf_d(I)+        &
     &      VWIND_FD(cf_c22(I,1),cf_c22(I,2))*cf_b(I)*cf_d(I) )
          RECORD_OUT(I,3) = cf_J(I)*(                                 &
     &      MSLP_FD(cf_c11(I,1),cf_c11(I,2))*cf_a(I)*cf_c(I)+         &
     &      MSLP_FD(cf_c21(I,1),cf_c21(I,2))*cf_b(I)*cf_c(I)+         &
     &      MSLP_FD(cf_c12(I,1),cf_c12(I,2))*cf_a(I)*cf_d(I)+         &
     &      MSLP_FD(cf_c22(I,1),cf_c22(I,2))*cf_b(I)*cf_d(I) )
        END DO
      if(myrank==0) then
      WRITE(16,*) 'UWIND_FD, min/max=', minval(UWIND_FD), maxval(UWIND_FD)
      WRITE(16,*) 'VWIND_FD, min/max=', minval(VWIND_FD), maxval(VWIND_FD)
      WRITE(16,*) 'MSLP_FD , min/max=', minval(MSLP_FD), maxval(MSLP_FD)
      WRITE(16,*) 'UWIND_FE, min/max=', minval(RECORD_OUT(:,1)), maxval(RECORD_OUT(:,1))
      WRITE(16,*) 'VWIND_FE, min/max=', minval(RECORD_OUT(:,2)), maxval(RECORD_OUT(:,2))
      WRITE(16,*) 'MSLP_FE,  min/max=', minval(RECORD_OUT(:,3)), maxval(RECORD_OUT(:,3))
      FLUSH(16)
      endif
      END SUBROUTINE READ_INTERP_NETCDF_CF
      
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_DIRECT(RECORD_IN,RECORD_OUT)
      USE NETCDF
      USE schism_glbl, only : npa, np_global,ipgl,rkind
      USE schism_glbl, only : cf_add_offset_uwind, cf_scale_factor_uwind
      USE schism_glbl, only : cf_add_offset_vwind, cf_scale_factor_vwind
      USE schism_glbl, only : cf_add_offset_pr, cf_scale_factor_pr
      USE schism_msgp, only : myrank

      IMPLICIT NONE
      INTEGER, INTENT(in)                      :: RECORD_IN
      REAL(rkind),dimension(npa,3),INTENT(out) :: RECORD_OUT
      REAL(rkind)			  :: TMPIN(np_global,3)
      INTEGER                             :: FID, ID, ISTAT, I, ITMP
      
      character (len = *), parameter :: CallFct="READ_NETCDF_DIRECT"
      character(len=100)       :: CHRERR

      ISTAT = NF90_OPEN('UVP_direct.nc', NF90_NOWRITE, FID)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)
      IF (ISTAT /= 0) THEN
      CHRERR = nf90_strerror(ISTAT)
      if(myrank==0)  WRITE(16,*) 'There is no UVP_direct.nc file!! ', TRIM(CHRERR)
      ENDIF

      ISTAT = NF90_inq_varid(FID, 'Uwind', ID)
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)
      IF (ISTAT /= 0) THEN
      CHRERR = nf90_strerror(ISTAT)
      if(myrank==0)  WRITE(16,*) 'There is no Uwind variable ', TRIM(CHRERR)
      TMPIN(:,1)=0
      ELSE
      ISTAT = NF90_GET_VAR(FID, ID, TMPIN(:,1), start = (/ 1, RECORD_IN /), count = (/ np_global, 1 /))
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)
      ENDIF

      ISTAT = NF90_inq_varid(FID, 'Vwind', ID)
      CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)
      IF (ISTAT /= 0) THEN
      CHRERR = nf90_strerror(ISTAT)
      if(myrank==0)  WRITE(16,*) 'There is no Vwind variable ', TRIM(CHRERR)
      TMPIN(:,2)=0
      ELSE
      ISTAT = NF90_GET_VAR(FID, ID, TMPIN(:,2), start = (/ 1, RECORD_IN /), count = (/ np_global, 1 /))
      CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)
      ENDIF

      ISTAT = NF90_inq_varid(FID, 'Pair', ID)
      CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)
      IF (ISTAT /= 0) THEN
      CHRERR = nf90_strerror(ISTAT)
      if(myrank==0)  WRITE(16,*) 'There is no Pair variable ', TRIM(CHRERR)
      TMPIN(:,3)=101325
      ELSE
      ISTAT = NF90_GET_VAR(FID, ID, TMPIN(:,3), start = (/ 1, RECORD_IN /), count = (/ np_global, 1 /))
      CALL GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)
      ENDIF
 
      ! MSLP in the WRF file is in mbar but have to use Pa -> multiply with 100.0
      do I=1, np_global
                if(ipgl(I)%rank==myrank) then
                itmp=ipgl(i)%id
                RECORD_OUT(itmp,1)=cf_add_offset_uwind + cf_scale_factor_uwind*TMPIN(I,1)
                RECORD_OUT(itmp,2)=cf_add_offset_vwind + cf_scale_factor_vwind*TMPIN(I,2)
                RECORD_OUT(itmp,3)=100*(cf_add_offset_pr+cf_scale_factor_pr*TMPIN(I,3))
                endif
      end do
      if(myrank==0) then
      WRITE(16,*) 'UWIND_FE, min/max = ', minval(RECORD_OUT(:,1)), maxval(RECORD_OUT(:,1))
      WRITE(16,*) 'VWIND_FE, min/max = ', minval(RECORD_OUT(:,2)), maxval(RECORD_OUT(:,2))
      WRITE(16,*) 'MSLP_FE,  min/max = ', minval(RECORD_OUT(:,3)), maxval(RECORD_OUT(:,3))
      FLUSH(16)
      endif

      ISTAT = NF90_CLOSE(FID)
      CALL GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)
      if(myrank==0) then
      WRITE(16,'(A, I6)') 'Done with UVP_direct.nc reading data for record = ', RECORD_IN
      FLUSH(16)
      endif
      END SUBROUTINE READ_NETCDF_DIRECT

!**********************************************************************
!*                                                                    *
!**********************************************************************      
      SUBROUTINE GENERIC_NETCDF_ERROR(CallFct, idx, iret)
      USE NETCDF

      implicit none
      integer, intent(in) :: iret, idx
      character(*), intent(in) :: CallFct
      character(len=500) :: CHRERR
      IF (iret .NE. nf90_noerr) THEN
        CHRERR = nf90_strerror(iret)
        WRITE(16,*) TRIM(CallFct), ' -', idx, '-', CHRERR
      ENDIF
      END SUBROUTINE GENERIC_NETCDF_ERROR

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CF_EXTRACT_TIME_NEW(eStrUnitTime, ConvertToSec)
      USE schism_glbl, only : rkind 
      USE schism_msgp, only : myrank
      IMPLICIT NONE
      character(len=100), intent(in) :: eStrUnitTime
      integer, intent(out) :: ConvertToSec
      character (len=100) :: Xname, Yname
      character (len=10) :: YnameYear, YnameMonth, YnameDay
      character (len=10) :: YnameHour, YnameMin, YnameSec
      character (len=50) :: YnameB, YnameD, YnameE
      character (len=50) :: YnameDate, YnameTime, YnameTimeP
      character (len=15) :: eStrTime
      integer alenB, alenC, alenD, alenE, alenTime, alenDate
      integer alen, posBlank
      integer lenHour, lenMin, lenSec, lenMonth, lenDay, posSepDateTime
      alen=LEN_TRIM(eStrUnitTime)
      posBlank=INDEX(eStrUnitTime(1:alen), ' ')
      Xname=eStrUnitTime(1:posBlank-1) ! should be days/hours/seconds
      IF (TRIM(Xname) .eq. 'days') THEN
        ConvertToSec=86400
      ELSEIF (TRIM(Xname) .eq. 'hours') THEN
        ConvertToSec=3600
      ELSEIF (TRIM(Xname) .eq. 'seconds') THEN
        ConvertToSec=1
      ELSE
        if(myrank==0)  WRITE(16,*) 'Error in the code for conversion'
      END IF
      !
      Yname=eStrUnitTime(posBlank+1:alen)
      alenB=LEN_TRIM(Yname)
      posBlank=INDEX(Yname(1:alenB), ' ')
      YnameB=Yname(posBlank+1:alenB) ! should be 1990-01-01 0:0:0
      !
      alenC=LEN_TRIM(YnameB)
      posSepDateTime=INDEX(YnameB(1:alenC), ' ')
      IF (posSepDateTime .gt. 0) THEN
        YnameDate=YnameB(1:posSepDateTime-1) ! should be 1990-01-01
        YnameTimeP=YnameB(posSepDateTime+1:alenC) ! should be 0:0:0
        alenC=LEN_TRIM(YnameTimeP)
        posBlank=INDEX(YnameTimeP(1:alenC), ' ')
        IF (posBlank .eq. 0) THEN
          YnameTime=YnameTimeP
        ELSE
          YnameTime=YnameTimeP(1:posBlank-1)
        END IF
      ELSE
        YnameDate=YnameB
        eStrTime(10:10)='0'
        eStrTime(11:11)='0'
        eStrTime(12:12)='0'
        eStrTime(13:13)='0'
        eStrTime(14:14)='0'
        eStrTime(15:15)='0'
      END IF
      !
      alenDate=LEN_TRIM(YnameDate)
      posBlank=INDEX(YnameDate(1:alenDate), '-')
      YnameYear=YnameDate(1:posBlank-1) ! should be 1990
      YnameD=YnameDate(posBlank+1:alenDate)
      alenD=LEN_TRIM(YnameD)
      posBlank=INDEX(YnameD(1:alenD), '-')
      YnameMonth=YnameD(1:posBlank-1) ! should be 01
      YnameDay=YnameD(posBlank+1:alenD) ! should be 01
      !
      ! year
      eStrTime( 1: 1)=YnameYear( 1: 1)
      eStrTime( 2: 2)=YnameYear( 2: 2)
      eStrTime( 3: 3)=YnameYear( 3: 3)
      eStrTime( 4: 4)=YnameYear( 4: 4)
      !
      ! month
      lenMonth=LEN_TRIM(YnameMonth)
      IF (lenMonth .eq. 2) THEN
        eStrTime( 5: 5)=YnameMonth( 1: 1)
        eStrTime( 6: 6)=YnameMonth( 2: 2)
      ELSE
        IF (lenMonth .eq. 1) THEN
          eStrTime( 5: 5)='0'
          eStrTime( 6: 6)=YnameMonth( 1: 1)
        ELSE
          if(myrank==0)  WRITE(16,*) 'DIE in trying to get the month'
        END IF
      END IF
      !
      ! day
      lenDay=LEN_TRIM(YnameDay)
      IF (lenDay .eq. 2) THEN
        eStrTime( 7: 7)=YnameDay( 1: 1)
        eStrTime( 8: 8)=YnameDay( 2: 2)
      ELSE
        IF (lenDay .eq. 1) THEN
          eStrTime( 7: 7)='0'
          eStrTime( 8: 8)=YnameDay( 1: 1)
        ELSE
          if(myrank==0)  WRITE(16,*) 'DIE in trying to get the day'
        END IF
      END IF
      !
      eStrTime( 9: 9)='.'
      !
      IF (posSepDateTime .gt. 0) THEN
        !
        alenTime=LEN_TRIM(YnameTime)
        posBlank=INDEX(YnameTime(1:alenTime), ':')
        YnameHour=YnameTime(1:posBlank-1) ! should be 0
        YnameE=YnameTime(posBlank+1:alenTime)
        alenE=LEN_TRIM(YnameE)
        posBlank=INDEX(YnameE(1:alenE), ':')
        YnameMin=YnameE(1:posBlank-1) ! should be 0
        YnameSec=YnameE(posBlank+1:alenE) ! should be 0
        !
        !
        ! Hour
        lenHour=LEN_TRIM(YnameHour)
        IF (lenHour .eq. 2) THEN
          eStrTime(10:10)=YnameHour( 1: 1)
          eStrTime(11:11)=YnameHour( 2: 2)
        ELSE
          IF (lenHour .eq. 1) THEN
            eStrTime(10:10)='0'
            eStrTime(11:11)=YnameHour( 1: 1)
          ELSE
            if(myrank==0)  WRITE(16,*) 'DIE in trying to get the hour'
          END IF
        END IF
        !
        ! Min
        lenMin=LEN_TRIM(YnameMin)
        IF (lenMin .eq. 2) THEN
          eStrTime(12:12)=YnameMin( 1: 1)
          eStrTime(13:13)=YnameMin( 2: 2)
        ELSE
          IF (lenMin .eq. 1) THEN
            eStrTime(12:12)='0'
            eStrTime(13:13)=YnameMin( 1: 1)
          ELSE
            if(myrank==0)  WRITE(16,*) 'DIE in trying to get the min'
          END IF
        END IF
        !
        ! Sec
        lenSec=LEN_TRIM(YnameSec)
        IF (lenSec .eq. 2) THEN
          eStrTime(14:14)=YnameSec( 1: 1)
          eStrTime(15:15)=YnameSec( 2: 2)
        ELSE
          IF (lenSec .eq. 1) THEN
            eStrTime(14:14)='0'
            eStrTime(15:15)=YnameSec( 1: 1)
          ELSE
            if(myrank==0)  WRITE(16,*) 'DIE in trying to get the sec'
          END IF
        END IF
      END IF
      if(myrank==0)  WRITE(16,*) 'WIND_START_TIME =', eStrTime
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_FRC_REC(wtime, record1, record2, w1, w2)
      ! For given WTIME return records to get and weights for time
      ! interpolation F(wwm_time)=F(rec1)*w1 + F(rec2)*w2
      ! wind_time is in the same units (s) as model_time 

      USE schism_glbl, only : nwtimes, wind_time_sec, rkind
      USE schism_msgp, only : myrank

      IMPLICIT NONE
      REAL(rkind), INTENT(IN)             :: wtime
      REAL(rkind), INTENT(OUT)            :: w1, w2
      INTEGER, INTENT(OUT)                :: record1, record2
      REAL(rkind) :: eTime1, eTime2
      INTEGER  :: iTime, I
 
      DO iTime=2,nwtimes
        eTime1=wind_time_sec(iTime-1)
        eTime2=wind_time_sec(iTime)
        IF ((eTime1 .le. wtime).and.(wtime .le. eTime2)) THEN
          record2=iTime
          record1=iTime-1
          w2=(wtime - eTime1)/(eTime2-eTime1)
          w1=(eTime2 - wtime)/(eTime2-eTime1)
          RETURN
        END IF
      END DO
	IF (record1.eq.0) THEN
	  record1=1
	  record2=2
	  w1=1
	  w2=0
	END IF
      END SUBROUTINE GET_FRC_REC
!**********************************************************************
!*                                                                    *
!**********************************************************************
