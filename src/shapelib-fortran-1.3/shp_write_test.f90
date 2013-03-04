PROGRAM shp_write_test
USE shplib
IMPLICIT NONE

INTEGER(kind=shp_ptr_c) :: shphandle
TYPE(shpobject),POINTER :: shpobj
INTEGER :: i ! , nent, tshp
CHARACTER(len=1024) :: filename

TYPE(shpobject),POINTER :: shpobj_new
REAL(kind=shp_fp_d) :: x(3),y(3),z(3)

CALL getarg(1,filename)
IF (filename == '') THEN
  PRINT'(A)','Usage: shape_test <shp_file>'
  STOP
ENDIF

shphandle = shpcreate(TRIM(filename), 5)

x(1)=0.6d0  ;  y(1)=0.3d0 ;  z(1)=0.d0
x(2)=0.9d0  ;  y(2)=0.3d0 ; z(2)=0.d0
x(3)=0.75d0  ;  y(3)=0.6d0 ; z(3)=0.d0

shpobj_new => shpcreatesimpleobject(5,3,x,y,z)
i=shpwriteobject(shphandle,-1,shpobj_new)

CALL shpclose(shphandle)

END PROGRAM shp_write_test

