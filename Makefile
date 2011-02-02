#-*- mode: makefile; mode: font-lock; vc-back-end: RCS -*-
#=============================================================================#
#                      H   S      A   L   K   A   N   E                       #
#=============================================================================#
#                                                                             #
# $Id: Makefile,v 1.1 2011/02/02 11:48:36 phseal Exp $
#                                                                             #
#-----------------------------------------------------------------------------#
# D. Quigley, University of Warwick                                           #
#
# $Log: Makefile,v $
# Revision 1.1  2011/02/02 11:48:36  phseal
# Initial revision
#
#
#
#=============================================================================#

SHELL = /bin/sh


# Where you want the executable
prefix     = $(HOME)
bindir     = $(prefix)/bin

# Compiler and flags
F90       = ifort
LD        = ifort
FFLAGS    = -O2 
INCLUDE   = 
LIBS      = 

# Nothing should need to be changed below here.
#=====================================================================
# Define objects in dependency order
OBJECTS   = constants.o timer.o random.o quaternion.o box.o vis_module.o alkane.o mc.o io.o main.o


# Main build target - builds the main executable
hs_alkane : $(OBJECTS)

	$(LD) -o $(bindir)/hs_alkane $(OBJECTS) $(FFLAGS) $(LIBS) 


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

