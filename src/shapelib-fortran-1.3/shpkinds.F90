#include "config.h"

MODULE shpkinds
IMPLICIT NONE
! Derived from Type_Kinds.f90 by: Paul van Delst, CIMSS/SSEC 12-Jun-2000
!                                 paul.vandelst <at> ssec.wisc.edu
! Copyright (C) 2000, 2004 Paul van Delst
! Modified by Davide Cesari, 2007, dcesari <at> arpa.emr.it
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

PRIVATE

INTEGER, PARAMETER, PUBLIC :: &
 shp_int_b    = SELECTED_INT_KIND(1), & ! Byte  integer
 shp_int_s    = SELECTED_INT_KIND(4), & ! Short integer
 shp_int_l    = SELECTED_INT_KIND(8)    ! Long  integer
INTEGER, PARAMETER :: &
 int_ll_t = SELECTED_INT_KIND(16)  ! LLong integer
INTEGER, PARAMETER, PUBLIC :: &
 shp_int_ll   = ( ( ( 1 + SIGN( 1, int_ll_t ) ) / 2 ) * int_ll_t ) + &
 ( ( ( 1 - SIGN( 1, int_ll_t ) ) / 2 ) * shp_int_l    )

INTEGER, PARAMETER, PUBLIC :: &
 shp_fp_s = SELECTED_REAL_KIND(6), & ! Single precision
 shp_fp_d = SELECTED_REAL_KIND(15)   ! Double precision
INTEGER, PARAMETER :: &
 fp_q_t = SELECTED_REAL_KIND(20) ! Quad precision
INTEGER, PARAMETER, PUBLIC  :: &
 shp_fp_q   = ( ( ( 1 + SIGN( 1, fp_q_t ) ) / 2 ) * fp_q_t ) + &
 ( ( ( 1 - SIGN( 1, fp_q_t ) ) / 2 ) * shp_fp_d )

INTEGER, PARAMETER, PUBLIC :: &
 shp_ptr_c = SIZEOF_PTR_C ! C pointer

END MODULE shpkinds
