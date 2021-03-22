# python used to launch the test suite.
# The packages netCDF4 and numpy are needed
export SCWWM_PY=/home/lmentaschi/usr/miniconda3/bin/python

# mpirun
export SCWWM_MPIRUN=mpirun.mpich

# schismWWM executables
export SCWWM_EXE=/home/lmentaschi/src/git/wwmv/wwm_menta_uost/src/schism.ex_VL_WWM
export SCWWM_COMBOUTPUT_EXE=/home/lmentaschi/src/git/wwmv/wwm_menta_uost/src/Utility/Combining_Scripts/combine_output11
export SCWWM_COMBHOTSTART_EXE=/home/lmentaschi/src/git/wwmv/wwm_menta_uost/src/Utility/Combining_Scripts/combine_hotstart7

# adding all the runtime paths in LD_LIBRARY_PATH (e.g. mpi, netcdf, hdf5 ...)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu/
