PROGRAM		= PLOT_results \
		  AltimeterComparison \
		  CREATE_sflux \
		  INTERPOL_field \
		  PLOT_diff_results

OBJECTS		= $(PROGRAM:%=%.o)

#CC              = clang++ -Weverything -std=c++11
CC		= g++ -std=c++11 -Wall -g -O2 

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CFLAGS 		= -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(NETCDF_CXX_ALLINC) $(GRIB_INCDIR)

LIBS		= $(NETCDF_CXX_LINK) $(GRIB_C_LINK)

%.o: %.cpp
	@echo
	@echo compile $<
	$(CC) $(OUTPUTFLAG) $@ $(CFLAGS) $(NOLINKFLAG) $<

all:	$(PROGRAM)

$(PROGRAM): 	$(OBJECTS)
		@echo
		@echo link $@
		$(CC) $(DEBUGFLAG) $@.o $(LIBS) -o $@

clean:;
	rm -f *.o
	rm -f $(PROGRAM)
