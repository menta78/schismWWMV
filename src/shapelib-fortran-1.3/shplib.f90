MODULE shplib_def
USE shpkinds
IMPLICIT NONE

! Fortran 90 bindings for shapelib (http://shapelib.maptools.org/)
! Copyright (C) Davide Cesari, 2007, dcesari <at> arpa.emr.it
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

INTEGER, PARAMETER :: &
 SHPT_NULL = 0, &
!  2D Shape Types (pre ArcView 3.x):
 shpt_point = 1, & ! Points
 shpt_arc = 3, & ! Arcs (Polylines, possible in parts)
 shpt_polygon = 5, & ! Polygons (possible in parts)
 shpt_multipoint = 8, & ! MultiPoint (related points)
!  3D Shape Types (may include "measure" values for vertices):
 shpt_pointz = 11, &
 shpt_arcz = 13, &
 shpt_polygonz = 15, &
 shpt_multipointz = 18, &
!  2D + Measure Types:
 shpt_pointm = 21, &
 shpt_arcm = 23, &
 shpt_polygonm = 25, &
 shpt_multipointm = 28, &
!  Complex (TIN-like) with Z, and Measure:
 shpt_multipatch = 31

TYPE shpobject
  INTEGER(kind=shp_ptr_c) :: shpobject_orig
  INTEGER :: nshptype, & ! Shape Type (SHPT_* - see list above)
   nshapeid, & ! Shape Number (-1 is unknown/unassigned)
   nparts ! # of Parts (0 implies single part with no info)
  INTEGER, POINTER :: panpartstart(:), & ! Start Vertex of part
   panparttype(:) ! Part Type (SHPP_RING if not SHPT_MULTIPATCH)
  INTEGER :: nvertices ! Vertex list 
  REAL(kind=shp_fp_d), POINTER ::  padfx(:), padfy(:), &
   padfz(:), padfm(:) ! (all zero if not provided)
  REAL(kind=shp_fp_d) :: & ! Bounds in X, Y, Z and M dimensions
   dfxmin, dfymin, dfzmin, dfmmin, dfxmax, dfymax, dfzmax, dfmmax
END TYPE shpobject

END MODULE shplib_def

MODULE shplib
USE shplib_def
IMPLICIT NONE

INTERFACE

  FUNCTION shpopen(pszShapeFile, pszAccess)
  USE shplib_def
  CHARACTER(len=*) :: pszShapeFile, pszAccess
  INTEGER(kind=shp_ptr_c) :: shpopen
  END FUNCTION shpopen

  SUBROUTINE shpgetinfo(hshp, entities, shapetype, minbound, maxbound)
  USE shplib_def
  INTEGER(kind=shp_ptr_c) :: hshp
  INTEGER :: entities, shapetype
  REAL(kind=shp_fp_d), OPTIONAL :: minbound(4), maxbound(4)
  END SUBROUTINE shpgetinfo

  FUNCTION shpreadobject_int(hshp, ishape, ftnobject)
  USE shplib_def
  INTEGER(kind=shp_ptr_c) :: hshp
  INTEGER :: ishape
  TYPE(shpobject) :: ftnobject
  INTEGER :: shpreadobject_int
  END FUNCTION shpreadobject_int

  SUBROUTINE shpclose(hshp)
  USE shplib_def
  INTEGER(kind=shp_ptr_c) :: hshp
  END SUBROUTINE shpclose

  FUNCTION shpcreate(pszShapeFile, nshapetype)
  USE shplib_def
  CHARACTER(len=*) :: pszShapeFile
  INTEGER :: nshapetype
  INTEGER(kind=shp_ptr_c) :: shpcreate
  END FUNCTION shpcreate

  FUNCTION shpcreatesimpleobject_int(nshptype, nvertices, padfx, padfy, padfz, &
   ftnobject)
  USE shplib_def
  INTEGER :: nshptype, nvertices
  REAL(kind=shp_fp_d) :: padfx(nvertices), padfy(nvertices)
  REAL(kind=shp_fp_d), OPTIONAL :: padfz(nvertices)
  TYPE(shpobject) :: ftnobject
  INTEGER :: shpcreatesimpleobject_int
  END FUNCTION shpcreatesimpleobject_int

  FUNCTION shpcreateobject_int(nshptype, ishape, nparts, &
   panpartstart, panparttype, &
   nvertices, padfx, padfy, padfz, padfm, ftnobject)
  USE shplib_def
  INTEGER :: nshptype, ishape, nparts, nvertices
  INTEGER :: panpartstart(nparts), panparttype(nparts)
  REAL(kind=shp_fp_d) :: padfx(nvertices), padfy(nvertices)
  REAL(kind=shp_fp_d), OPTIONAL :: padfz(nvertices), padfm(nvertices)
  TYPE(shpobject) :: ftnobject
  INTEGER :: shpcreateobject_int
  END FUNCTION shpcreateobject_int

  SUBROUTINE shpcomputeextents_int(psobject, ftnobject)
  USE shplib_def
  INTEGER(kind=shp_ptr_c) :: psobject
  TYPE(shpobject) :: ftnobject
  END SUBROUTINE shpcomputeextents_int
  
  FUNCTION shpwriteobject_int(hshp, ishape, psobject)
  USE shplib_def
  INTEGER(kind=shp_ptr_c) :: hshp
  INTEGER :: ishape
  INTEGER(kind=shp_ptr_c) :: psobject
  INTEGER :: shpwriteobject_int
  END FUNCTION shpwriteobject_int

  SUBROUTINE shpdestroyobject_int(psobject)
  USE shplib_def
  INTEGER(kind=shp_ptr_c) :: psobject
  END SUBROUTINE shpdestroyobject_int

  FUNCTION shprewindobject_int(hshp, psobject)
  USE shplib_def
  INTEGER(kind=shp_ptr_c) :: hshp
  INTEGER(kind=shp_ptr_c) :: psobject
  INTEGER :: shprewindobject_int
  END FUNCTION shprewindobject_int

END INTERFACE

CONTAINS


FUNCTION shpreadobject(hshp, ishape)
INTEGER(kind=shp_ptr_c) :: hshp
INTEGER :: ishape
TYPE(shpobject), POINTER :: shpreadobject

INTEGER :: ier

ALLOCATE(shpreadobject)
ier = shpreadobject_int(hshp, ishape, shpreadobject)
IF (ier /= 0) DEALLOCATE(shpreadobject)

END FUNCTION shpreadobject


FUNCTION shpcreatesimpleobject(nshptype, nvertices, padfx, padfy, padfz)
INTEGER :: nshptype, nvertices
REAL(kind=shp_fp_d) :: padfx(nvertices), padfy(nvertices)
REAL(kind=shp_fp_d), OPTIONAL :: padfz(nvertices)
TYPE(shpobject), POINTER :: shpcreatesimpleobject

INTEGER :: ier

ALLOCATE(shpcreatesimpleobject)
ier = shpcreatesimpleobject_int(nshptype, nvertices, padfx, padfy, padfz, &
 shpcreatesimpleobject)
IF (ier /= 0) DEALLOCATE(shpcreatesimpleobject)

END FUNCTION shpcreatesimpleobject


FUNCTION shpcreateobject(nshptype, ishape, nparts, panpartstart, panparttype, &
 nvertices, padfx, padfy, padfz, padfm)
INTEGER :: nshptype, ishape, nparts, nvertices
INTEGER :: panpartstart(nparts), panparttype(nparts)
REAL(kind=shp_fp_d) :: padfx(nvertices), padfy(nvertices)
REAL(kind=shp_fp_d), OPTIONAL :: padfz(nvertices), padfm(nvertices)
TYPE(shpobject), POINTER :: shpcreateobject

INTEGER :: ier

ALLOCATE(shpcreateobject)
ier = shpcreateobject_int(nshptype, ishape, nparts, panpartstart, panparttype, &
 nvertices, padfx, padfy, padfz, padfm, shpcreateobject)
IF (ier /= 0) DEALLOCATE(shpcreateobject)

END FUNCTION shpcreateobject


SUBROUTINE shpcomputeextents(psobject)
TYPE(shpobject) :: psobject

CALL shpcomputeextents_int(psobject%shpobject_orig, psobject)

END SUBROUTINE shpcomputeextents


FUNCTION shpwriteobject(hshp, ishape, psobject)
INTEGER(kind=shp_ptr_c) :: hshp
INTEGER :: ishape
TYPE(shpobject) :: psobject
INTEGER :: shpwriteobject

shpwriteobject = shpwriteobject_int(hshp, ishape, psobject%shpobject_orig)

END FUNCTION shpwriteobject


SUBROUTINE shpdestroyobject(psobject)
TYPE(shpobject), POINTER :: psobject

IF (ASSOCIATED(psobject)) THEN
  CALL shpdestroyobject_int(psobject%shpobject_orig)
!?!?!  CALL shpcomputeextents(psobject)
  DEALLOCATE(psobject)
ENDIF

END SUBROUTINE shpdestroyobject


FUNCTION shprewindobject(hshp, psobject)
INTEGER(kind=shp_ptr_c) :: hshp
TYPE(shpobject) :: psobject
LOGICAL :: shprewindobject

INTEGER :: ier

ier = shprewindobject_int(hshp, psobject%shpobject_orig)
IF (ier == 0) THEN
  shprewindobject = .FALSE.
ELSE
  shprewindobject = .TRUE.
ENDIF
END FUNCTION shprewindobject


END MODULE shplib


SUBROUTINE shpset_object(obj, obj_orig, nshptype, nshapeid, &
 nparts, panpartstart, panparttype, &
 nvertices, padfx, padfy, padfz, padfm, &
 dfxmin, dfymin, dfzmin, dfmmin, dfxmax, dfymax, dfzmax, dfmmax)
USE shplib_def
IMPLICIT NONE

TYPE(shpobject) :: obj
INTEGER(kind=shp_ptr_c) :: obj_orig
INTEGER :: nshptype ! Shape Type (SHPT_* - see list above)
INTEGER :: nshapeid ! Shape Number (-1 is unknown/unassigned)
INTEGER :: nparts ! # of Parts (0 implies single part with no info)
INTEGER, TARGET :: panpartstart(nparts), & ! Start Vertex of part
 panparttype(nparts) ! Part Type (SHPP_RING if not SHPT_MULTIPATCH)
INTEGER :: nvertices ! Vertex list 
REAL(kind=shp_fp_d), TARGET ::  padfx(nvertices), padfy(nvertices), &
 padfz(nvertices), padfm(nvertices) ! (all zero if not provided)
REAL(kind=shp_fp_d) :: & ! Bounds in X, Y, Z and M dimensions
 dfxmin, dfymin, dfzmin, dfmmin, dfxmax, dfymax, dfzmax, dfmmax

obj%shpobject_orig = obj_orig
obj%nshptype = nshptype
obj%nshapeid = nshapeid
obj%nparts = nparts
obj%panpartstart => panpartstart
obj%panparttype => panparttype
obj%nvertices = nvertices
obj%padfx => padfx
obj%padfy => padfy
obj%padfz => padfz
obj%padfm => padfm
obj%dfxmin = dfxmin
obj%dfymin = dfymin
obj%dfzmin = dfzmin
obj%dfmmin = dfmmin
obj%dfxmax = dfxmax
obj%dfymax = dfymax
obj%dfzmax = dfzmax
obj%dfmmax = dfmmax

END SUBROUTINE shpset_object
