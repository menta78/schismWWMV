PROGRAM shp_read_test
USE shplib
IMPLICIT NONE

INTEGER(kind=shp_ptr_c) :: shphandle
REAL(kind=shp_fp_d) :: minbound(4), maxbound(4)
TYPE(shpobject),POINTER :: shpobj
INTEGER :: i, nent, tshp
CHARACTER(len=1024) :: filename

CALL getarg(1,filename)
IF (filename == '') THEN
  PRINT'(A)','Usage: shape_test <shp_file>'
  STOP
ENDIF
shphandle = shpopen(TRIM(filename), 'rb') ! TRIM(filename) is required
IF (shphandle == 0) THEN
  PRINT'(2A)','Error, cannot open shape file ',TRIM(filename)
  STOP
ENDIF

CALL shpgetinfo(shphandle, nent, tshp, minbound, maxbound)
PRINT'(2(A,I0))','Number of shapes: ',nent,' Type of shapes: ',tshp
DO i = 0, nent-1
  shpobj => shpreadobject(shphandle, i)
  CALL shpcomputeextents(shpobj)
  CALL shp_print(shpobj)
  CALL shpdestroyobject(shpobj)
ENDDO

CALL shpclose(shphandle)

CONTAINS

SUBROUTINE shp_print(shpobj)
TYPE(shpobject) :: shpobj

INTEGER :: i

PRINT'(4(A,I0))','Shape ',shpobj%nshapeid,'  Type: ',shpobj%nshptype, &
 ' No. of parts: ',shpobj%nparts,' No. of vertices: ',shpobj%nvertices
IF (shpobj%nvertices /= SIZE(shpobj%padfx)) THEN
  PRINT*,shpobj%nvertices,SIZE(shpobj%padfx)
  RETURN
ENDIF
PRINT'(A,4G12.6)','Shape limits: ',shpobj%dfxmin,shpobj%dfymin,shpobj%dfxmax,&
 shpobj%dfymax
DO i = 1, shpobj%nvertices
  PRINT'(2G12.6)',shpobj%padfx(i),shpobj%padfy(i)
ENDDO

END SUBROUTINE shp_print

END PROGRAM shp_read_test

