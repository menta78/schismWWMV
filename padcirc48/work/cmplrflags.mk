# SRCDIR is set in makefile or on the compile line
INCDIRS := -I . -I $(SRCDIR)/prep
########################################################################
# Compiler flags for Linux operating system on x86_64-unknown-linux-gnu
# IVICA BOREAS x86_64-unknown-linux-gnu
#
compiler=intel
# Intel compiler
ifeq ($(compiler),intel) 
  PPFC	        :=  ifort -w  
  FC	        :=  ifort -w   
  PFC	        :=  mpif90
  OPTLVL        := -O2 
  ifeq ($(ADC_DEBUG),yes)
    OPTLVL        :=  
  endif
  FFLAGS1	:=  $(INCDIRS) $(OPTLVL) -extend_source -assume byterecl -g -check bounds -traceback  
  FFLAGS2	:=  $(FFLAGS1)
  FFLAGS3	:=  $(FFLAGS1)
  DA  	        :=  -DREAL8 -DLINUX -DCSCA
  DP  	        :=  -DREAL8 -DLINUX -DCSCA -DCMPI -DHAVE_MPI_MOD
  DPRE	        :=  -DREAL8 -DLINUX
  IMODS 	:=  -I
  CC            := mpicc
  CCBE          := mpicc	
  CFLAGS        := $(INCDIRS) $(OPTLVL) -DLINUX
  CLIBS         := 
  LIBS  	:=  
  MSGLIBS	:=  
  $(warning (INFO) Corresponding machine found in cmplrflags.mk.)
  ifneq ($(FOUND),TRUE)
     FOUND := TRUE
  else 
     MULTIPLE := TRUE
  endif
endif
