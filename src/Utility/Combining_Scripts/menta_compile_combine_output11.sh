module purge

module load gcc/12.2.0
module load netcdf-fortran/4.6.1--gcc--12.2.0 # should also load te dependent modules

gfortran -O3 -Bstatic -o combine_output11 -ffree-line-length-none ../UtilLib/schism_geometry.f90 ../UtilLib/argparse.f90 combine_output11.f90 -I/leonardo/prod/spack/5.2/install/0.21/linux-rhel8-icelake/gcc-12.2.0/netcdf-fortran-4.6.1-ioyufgaugvqjspshyaje65wscwubpvpz/include/ -I/leonardo/prod/spack/5.2/install/0.21/linux-rhel8-icelake/gcc-12.2.0/netcdf-c-4.9.2-wk7h7wank6kyxvcq7gru5taerb2dec3t/include -lnetcdff -lnetcdf

