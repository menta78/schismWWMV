#!/bin/bash
callpth=${BASH_SOURCE[0]}
callpth=$(readlink -f $callpth)
bashsrcdir=$(cd `dirname "$callpth"` && pwd)

module load netcdf-fortran/4.6.1--gcc--12.2.0 # should also load te dependent modules

$bashsrcdir/combine_hotstart7 $@

