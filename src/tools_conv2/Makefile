PROG1 = $(HOME)/bin/ConvertBin_to_NC

OBJS1 = ConvertBin_to_NC.o
LIBS = -L$(NETCDF_LIBDIR) -lnetcdf

F90 = ifort

F90FLAGS = -traceback -check uninit  -warn interfaces,nouncalled -gen-interfaces -O1

NETC = NETC

$(NETC):
	cp -f $(NETCDF_INCDIR)/netcdf.mod .
	cp -f $(NETCDF_INCDIR)/typesizes.mod .

all: $(NETC) $(PROG1)

$(PROG1): $(OBJS1)
	$(F90) -o $(PROG1) $(F90FLAGS) $(OBJS1) $(LIBS) $(LDFLAGS)

clean:
	rm -f *.o *.mod *_mod.* *_genmod.* $(PROG)

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS)  -c $<
