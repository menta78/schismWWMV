PROGRAM = AltimeterComparison


OBJECTS		= $(PROGRAM:%=%.o)

#CC		= icpc -std=c++11 -Wall -g -O2 
#CC             = clang++ -Weverything -std=c++11
CC		= g++ -std=c++11 -Wall -g -O2 

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CFLAGS 		= -I$(EIGEN_PATH) $(NETCDF_CXX_ALLINC) $(GRIB_INCDIR)

LIBS		= -lm -lpthread $(NETCDF_CXX_LINK) $(GRIB_C_LINK)

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
