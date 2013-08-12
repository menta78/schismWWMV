# Copyright 2005-2013 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
#
# In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
# virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
set -ea

echo
echo "TEST: $0"

if [ -z "${data_dir}" ]
then
  cd ../../
  cpath=`pwd`
  GRIB_DEFINITION_PATH=$cpath/definitions
  export GRIB_DEFINITION_PATH
  GRIB_SAMPLES_PATH=$cpath/samples
  export GRIB_SAMPLES_PATH
  tools_dir=$cpath/tools/
  examples_dir=$cpath/examples/F90/
  data_dir=$cpath/data
  samples_dir=$cpath/samples
else
  echo "Skipping test $0"
  exit
fi

cd ${examples_dir}

if [ -z "${GRIB_API_INCLUDE}" ]
then 
  GRIB_API_INCLUDE=`pwd`/src
fi

if [ -z "${GRIB_API_LIB}" ]
then 
  GRIB_API_LIB=`pwd`/src
fi

#${tools_dir}grib_info

set -u

