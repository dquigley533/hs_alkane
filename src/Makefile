#=============================================================================#
#                      H   S      A   L   K   A   N   E                       #
#=============================================================================#
#                                                                             #
# $Id: Makefile,v 1.8 2012/07/20 12:40:00 phrkao Exp $
#                                                                             #
#-----------------------------------------------------------------------------#
# D. Quigley, University of Warwick                                           #
#=============================================================================#
SHELL = /bin/sh

# Where you want the executable
prefix     = $(HOME)
bindir     = $(prefix)/bin
libdir     = $(prefix)/lib

# Default compiler and flags
F90       = gfortran
LD        = gfortran
FFLAGS    = -O3 -fPIC -g 
INCLUDE   = 
LIBS      = 

# Nothing should need to be changed below here.
#=====================================================================
# Define objects in dependency order for stand-alone fortran90 code
OBJECTS   = constants.o timer.o random.o quaternion.o box.o alkane.o vis_module.o mc.o io.o main.o

# Define objects in dependency order for library backend to other codes
LIBOBJ = constants.o timer.o random.o quaternion.o box.o alkane.o vis_module.o mc_dummy.o io.o 

# Main build target - builds the main executable
hs_alkane : $(OBJECTS)

	$(LD) -o hs_alkane $(OBJECTS) $(FFLAGS) $(LIBS) 
	cp hs_alkane $(bindir)/hs_alkane

# Library - should make clean && make if building this after the stand-alone f90 code
library: $(LIBOBJ)
	$(F90) $(LIBOBJ) -shared -o $(libdir)/libalkane.so $(FFLAGS) $(LIBS) 

# Lists the procedues in libalkane.so accessible by C binding without mangling
library-query: library
	nm $(libdir)/libalkane.so | grep -v '__' | grep T

# Runs the single chain example
singlechain: hs_alkane
	$(bindir)/hs_alkane ../examples/singlechain.input

# Runs the solid butane example
solid_butane: hs_alkane 
	ln -s ../examples/chain.xmol .
	$(bindir)/hs_alkane ../examples/solid_butane.input

# Runs the liquid butane example
liquid_butane: hs_alkane
	$(bindir)/hs_alkane ../examples/liquid_butane.input

.PRECIOUS: %.o
.PHONY:  clean

%: %.o
%.o: %.f90
	$(F90) $(FFLAGS) $(INCLUDE) -c -o $@ $<

%.o: %.F90
	$(F90) $(FFLAGS) $(INCLUDE) -c -o $@ $<

%.o: %.f95
	$(F90) $(FFLAGS) $(INCLUDE) -c -o $@ $<
%: %.o
	$(F90) $(FFLAGS) $(INCLUDE) -o $@ $^

clean : 

	rm -f *.mod *.d *.il *.o work.*
	rm -f $(bindir)/hs_alkane $(libdir)/libalkane.so

