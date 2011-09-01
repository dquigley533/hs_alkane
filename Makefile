#-*- mode: makefile; mode: font-lock; vc-back-end: RCS -*-
#=============================================================================#
#                      H   S      A   L   K   A   N   E                       #
#=============================================================================#
#                                                                             #
# $Id: Makefile,v 1.5 2011/09/01 16:55:13 phrkao Exp $
#                                                                             #
#-----------------------------------------------------------------------------#
# D. Quigley, University of Warwick                                           #
#
# $Log: Makefile,v $
# Revision 1.5  2011/09/01 16:55:13  phrkao
# Changed alkane_check_chain_geometry to be C compatible, had to change
# the argument "violate" from logical to integer and subsequently changed
# mc.f90 where this was used.
#
# Revision 1.4  2011/08/02 12:56:47  phseal
# Added C bindings to all procedures which should be callable externally
# when compiled as a library.
#
# Revision 1.3  2011/08/02 10:55:11  phseal
# Initial version for compilation as a library
#
# Revision 1.2  2011/07/29 15:58:29  phseal
# Added multiple simulation box support.
#
# Revision 1.1.1.1  2011/02/02 11:48:36  phseal
# Initial import from prototype code.
#
#
#
#=============================================================================#

SHELL = /bin/sh


# Where you want the executable
prefix     = $(HOME)
bindir     = $(prefix)/bin

# Compiler and flags
F90       = gfortran-mp-4.5
LD        = gfortran-mp-4.5
FFLAGS    = -O0 -g -fPIC
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

	$(LD) -o $(bindir)/hs_alkane $(OBJECTS) $(FFLAGS) $(LIBS) 

# Library - should make clean && make if building this after the stand-along f90 code
library: $(LIBOBJ)
	$(F90) $(LIBOBJ) -shared -o libalkane.so $(FFLAGS) $(LIBS) 

# Lists the procedues in libalkane.so accessible by C binding without mangling
library-query:
	nm libalkane.so | grep -v '__' | grep T


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
	rm -f $(bindir)/hs_alkane

