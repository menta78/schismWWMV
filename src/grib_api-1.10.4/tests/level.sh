#!/bin/sh
# Copyright 2005-2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
# virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
#


. ./include.sh
#set -x

outfile=out.grib

file=${data_dir}/regular_gaussian_pressure_level.grib1

cat >good<<EOF
grib1: level=850 isobaricInhPa indicatorOfTypeOfLevel=100
grib2: level=850 isobaricInhPa pl 85000
grib1: level=850 isobaricInhPa indicatorOfTypeOfLevel=100
grib1: level=50 isobaricInPa indicatorOfTypeOfLevel=210
grib2: level=50 isobaricInPa pl 50
grib1: level=50 isobaricInPa indicatorOfTypeOfLevel=210
grib1: level=1 topLevel=1 bottomLevel=2 depthBelowLandLayer indicatorOfTypeOfLevel=112
grib2: level=1 topLevel=1 bottomLevel=2 depthBelowLandLayer first= sfc 1 second= sfc 2
grib1: level=1 topLevel=1 bottomLevel=2 depthBelowLandLayer indicatorOfTypeOfLevel=112
EOF

cat >level.filter<<EOF
set typeOfLevel="isobaricInhPa";
set level=850;
print "grib[edition]: level=[level] [typeOfLevel] indicatorOfTypeOfLevel=[indicatorOfTypeOfLevel:l]";
set edition=2;
print "grib[edition]: level=[level] [typeOfLevel] [typeOfFirstFixedSurface:s] [scaledValueOfFirstFixedSurface]";
set edition=1;
print "grib[edition]: level=[level] [typeOfLevel] indicatorOfTypeOfLevel=[indicatorOfTypeOfLevel:l]";
set edition=1;
set typeOfLevel="isobaricInPa";
set level=50;
print "grib[edition]: level=[level] [typeOfLevel] indicatorOfTypeOfLevel=[indicatorOfTypeOfLevel:l]";
set edition=2;
print "grib[edition]: level=[level] [typeOfLevel] [typeOfFirstFixedSurface:s] [scaledValueOfFirstFixedSurface]";
set edition=1;
print "grib[edition]: level=[level] [typeOfLevel] indicatorOfTypeOfLevel=[indicatorOfTypeOfLevel:l]";
set typeOfLevel="depthBelowLandLayer";
set topLevel=1;
set bottomLevel=2;
print "grib[edition]: level=[level] topLevel=[topLevel] bottomLevel=[bottomLevel] [typeOfLevel] indicatorOfTypeOfLevel=[indicatorOfTypeOfLevel:l]";
set edition=2;
print "grib[edition]: level=[level] topLevel=[topLevel] bottomLevel=[bottomLevel] [typeOfLevel] first= [typeOfFirstFixedSurface:s] [scaledValueOfFirstFixedSurface] second= [typeOfSecondFixedSurface:s] [scaledValueOfSecondFixedSurface]";
set edition=1;
print "grib[edition]: level=[level] topLevel=[topLevel] bottomLevel=[bottomLevel] [typeOfLevel] indicatorOfTypeOfLevel=[indicatorOfTypeOfLevel:l]";
EOF

    
${tools_dir}grib_filter level.filter $file > test.dump
diff good test.dump

rm -f level.filter good test.dump


