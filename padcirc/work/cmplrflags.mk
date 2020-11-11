# SRCDIR is set in makefile or on the compile line
INCDIRS := -I . -I/home/aron/padcirc/prep 
DEBUG=full
########################################################################
# Compiler flags for Linux operating system on 64bit x86 CPU
#
ifeq ($(MACHINE)-$(OS),x86_64-linux-gnu)
#
# ***NOTE*** User must select between various Linux setups
#            by commenting/uncommenting the appropriate compiler
#
compiler=ncep
#
# Compiler Flags for gfortran and gcc
ifeq ($(compiler),ncep)
  PPFC          := mpif90#${COMP} 
  FC            := ifort#${COMP}
  PFC           := mpif90#${COMP_MPI} 
  INCDIRS       := $(INCDIRS) -I${NETCDF_INCDIR} ${HDF5_INCDIR} -I${Z_INC}
  FFLAGS1       :=  $(INCDIRS) -FI -assume byterecl -132 -assume buffered_io -fp-model strict -traceback -check all -warn interfaces,nouncalled -gen-interface
  ifeq ($(DEBUG),full)
     FFLAGS1       :=  $(INCDIRS) -g -O0 -traceback -check all -warn interfaces,nouncalled -gen-interface -FI -assume byterecl -132 -DALL_TRACE -DFULL_STACK -DFLUSH_MESSAGES
  endif
  FFLAGS2       :=  $(FFLAGS1)
  FFLAGS3       :=  $(FFLAGS1)
  DA            :=  -DREAL8 -DLINUX -DCSCA
  DP            :=  -DREAL8 -DLINUX -DCSCA -DCMPI
  DPRE          :=  -DREAL8 -DLINUX -DADCSWAN
  ifeq ($(SWAN),enable)
     DPRE          := $(DPRE) -DADCSWAN
  endif
  IMODS         :=  -I
  CC            := icc#${C_COMP} 
  CCBE		      := mpicc#${C_COMP_MP} 
  CFLAGS        := $(INCDIRS) -m64 -DLINUX -fp-model strict
  ifeq ($(DEBUG),full)
     CFLAGS        := $(INCDIRS) -g -O0 -march=k8 -m64 -mcmodel=medium -DLINUX
  endif
#  CLIBS         :=
#  FLIBS         :=
#  MSGLIBS       :=
#  NETCDFHOME    :=/usrx/local/prod/NetCDF/4.2/intel/haswell
#  HDF5HOME      :=/usrx/local/prod/HDF5/1.8.9/serial/intel/haswell/lib
#  ZHOME         :=/usrx/local/prod/zlib/1.2.7/intel/haswell/lib
  ifeq ($(NETCDF),enable)
     FLIBS          := $(FLIBS) -L${NETCDF_LIBDIR} -L${NETCDF_LIBDIR} ${HDF5_LDFLAGS} ${Z_LIB}
  endif
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else
     MULTIPLE := TRUE
  endif
endif
#
endif

