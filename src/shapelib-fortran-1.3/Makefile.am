lib_LTLIBRARIES = libshp-fortran.la

noinst_PROGRAMS = shp_write_test shp_read_test

libshp_fortran_la_SOURCES = shpkinds.F90 shplib.f90 shplibf.c
nodist_include_HEADERS = shpkinds.mod shplib_def.mod shplib.mod

shplib_def.$(OBJEXT): shpkinds.$(OBJEXT)
shplib_def.lo: shpkinds.$(OBJEXT)

shp_read_test_SOURCES = shp_read_test.f90
shp_write_test_SOURCES = shp_write_test.f90
LDADD = libshp-fortran.la

dist_pkgdoc_DATA = LICENSE shapelib-fortran.html

EXTRA_DIST = shapelib-fortran.spec

mostlyclean-compile:
	-rm -f *.$(OBJEXT) *.mod
