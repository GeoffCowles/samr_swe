##
## File:	$URL: file:///usr/casc/samrai/repository/SAMRAI/tags/v-2-2-1/examples/swe/Makefile.in $
## Package:	SAMRAI applications
## Copyright:	(c) 1997-2007 Lawrence Livermore National Security, LLC
## Revision:	$LastChangedRevision: 1704 $
## Modified:	$LastChangedDate: 2007-11-13 16:32:40 -0800 (Tue, 13 Nov 2007) $
## Description:	makefile for swe gas dynamics sample application
##

#SAMRAI        = ../../../../../SAMRAI
#SRCDIR        = ../../../../../SAMRAI/source/test/applications/Euler
#SUBDIR        = source/test/applications/Euler
#VPATH         = ../../../../../SAMRAI/source/test/applications/Euler
#TESTTOOLS     = ../../testtools
#OBJECT        = /opt/samrai_3.1.0/ifortgcc64_opt
#SAMRAI        = /opt/samrai_3.1.0/ifortgcc64_opt
OBJECT        = /opt/samrai_3.1.0/ifortgcc64_debug
SAMRAI        = /opt/samrai_3.1.0/ifortgcc64_debug

SRCDIR        = ./swe

default:        swe

#include Makefile.config
include $(SAMRAI)/config/Makefile.config

CPPFLAGS_EXTRA= -DDISPLAY -DNDIM=2  -DTESTING=0

CXX_OBJS      = main.o swe.o 
F_OBJS        =  gparms.o cntrl.o philim.o limiter.o rpn.o rpt.o flux2fw_geo.o  flux.o sedflux.o friction.o drycheck.o wd.o grad.o init.o setbathy.o calcdt.o c2f.o junkprobe.o

main:
		if test -f stamp-3d; then $(MAKE) clean; fi
		touch stamp-2d
		$(MAKE) PDIM=2 swe

swe:	$(CXX_OBJS) $(F_OBJS) $(LIBSAMRAIDEPEND)
		$(CXX) $(CXXFLAGS) $(LDFLAGS) $(CXX_OBJS) $(F_OBJS)	\
		$(LIBSAMRAI2D) $(LIBSAMRAI) $(LDLIBS) -o swe 

clean:
		$(RM) *.f swe *.o *.mod

redo:
		$(RM) swe core *.ii *.int.c
		$(RM) -r ti_files ii_files

include Makefile.depend

FORTRAN       = fortran

gparms.o:	$(FORTRAN)/gparms.f90
		$(F90) $(FFLAGS) -c $(FORTRAN)/gparms.f90 -o $@

cntrl.o:	$(FORTRAN)/cntrl.f90
		$(F90) $(FFLAGS) -c $(FORTRAN)/cntrl.f90 -o $@
		
philim.o:	$(FORTRAN)/philim.f
	    $(F90) $(FFLAGS) -c $(FORTRAN)/philim.f -o $@

limiter.o:	$(FORTRAN)/limiter.f
		$(F77) $(FFLAGS) -c $(FORTRAN)/limiter.f -o $@
		
refine.o:		$(FORTRAN)/refine.f90
		$(F90) $(FFLAGS) -c $(FORTRAN)/refine.f90 -o $@

flux.o:		$(FORTRAN)/flux.f90
		$(F90) $(FFLAGS) -c $(FORTRAN)/flux.f90 -o $@
		
flux2fw_geo.o:		$(FORTRAN)/flux2fw_geo.f
			$(F77) $(FFLAGS) -c $(FORTRAN)/flux2fw_geo.f -o $@

rpn.o:		$(FORTRAN)/rpn.f
		$(F77) $(FFLAGS) -c $(FORTRAN)/rpn.f -o $@
		
rpt.o:		$(FORTRAN)/rpt.f
		$(F77) $(FFLAGS) -c $(FORTRAN)/rpt.f -o $@

wd.o:		$(FORTRAN)/wd.f90
		$(F90) $(FFLAGS) -c $(FORTRAN)/wd.f90 -o $@

friction.o:	$(FORTRAN)/friction.f90
		$(F90) $(FFLAGS) -c $(FORTRAN)/friction.f90 -o $@
		
sedflux.o:	$(FORTRAN)/sedflux.f90
		$(F90) $(FFLAGS) -c $(FORTRAN)/sedflux.f90 -o $@
		
drycheck.o:	$(FORTRAN)/drycheck.f90
		$(F90) $(FFLAGS) -c $(FORTRAN)/drycheck.f90 -o $@
		
junkprobe.o: $(FORTRAN)/junkprobe.f90
		$(F90) $(FFLAGS) -c $(FORTRAN)/junkprobe.f90 -o $@

grad.o:		$(FORTRAN)/grad.f90
		$(F90) $(FFLAGS) -c $(FORTRAN)/grad.f90 -o $@

init.o:		$(FORTRAN)/init.f90
		$(F90) $(FFLAGS) -c $(FORTRAN)/init.f90 -o $@
		
setbathy.o:		$(FORTRAN)/setbathy.f90
		$(F90) $(FFLAGS) -c $(FORTRAN)/setbathy.f90 -o $@

calcdt.o:	$(FORTRAN)/calcdt.f90
		$(F90) $(FFLAGS) -c $(FORTRAN)/calcdt.f90 -o $@

c2f.o:		$(FORTRAN)/c2f.f90
		$(F90) $(FFLAGS) -c $(FORTRAN)/c2f.f90 -o $@

