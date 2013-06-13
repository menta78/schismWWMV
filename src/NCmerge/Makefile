PROG1 = $(HOME)/bin/MergeWWMhistory_NC
PROG2 = $(HOME)/bin/MergeWWMhotfile_NC
PROG3 = $(HOME)/bin/MergeWWMstation_NC

OBJScommon = wwm_datapl.o wwm_netcdf.o wwm_gridcf.o wwm_mjdv2.o wwm_aux.o wwm_input.o
OBJhistory = MergeWWMhistory_NC.o
OBJhotfile = MergeWWMhotfile_NC.o
OBJstation = MergeWWMstation_NC.o
LIBS = -L$(NETCDF_LIBDIR) -lnetcdf

F90 = ifort
#F90 = gfortran

F90FLAGS = -traceback -check uninit  -warn interfaces,nouncalled -gen-interfaces -O1 -DNCDF -DMERGE_OPERATION
#F90FLAGS = -g -O3 -ffree-line-length-0 -fbacktrace -Wall

NETC = NETC

$(NETC):
	cp -f $(NETCDF_INCDIR)/netcdf.mod .
	cp -f $(NETCDF_INCDIR)/typesizes.mod .

all: $(NETC) $(PROG1) $(PROG2) $(PROG3)

$(PROG1): $(OBJScommon) $(OBJhistory)
	$(F90) -o $(PROG1) $(F90FLAGS) $(OBJScommon) $(OBJhistory) $(LIBS) $(LDFLAGS)

$(PROG2): $(OBJScommon) $(OBJhotfile)
	$(F90) -o $(PROG2) $(F90FLAGS) $(OBJScommon) $(OBJhotfile) $(LIBS) $(LDFLAGS)

$(PROG3): $(OBJScommon) $(OBJstation)
	$(F90) -o $(PROG3) $(F90FLAGS) $(OBJScommon) $(OBJstation) $(LIBS) $(LDFLAGS)

clean:
	rm -f *.o *.mod *_mod.* *_genmod.* $(PROG)

.SUFFIXES: $(SUFFIXES) .F90

.F90.o:
	$(F90) $(F90FLAGS)  -c $<
