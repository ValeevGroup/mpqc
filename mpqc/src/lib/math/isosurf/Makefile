TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

LD=$(CXX)

GENCXXSRC = edgeSet.cc edgeAVLSet.cc \
            edgeMap.cc edgeRAVLMap.cc
GENINC = $(GENCXXSRC:%.cc=%.h)
CXXSRC = volume.cc shape.cc surf.cc volume_i.cc \
         vertex.cc edge.cc triangle.cc isosurf.cc
CSRC = implicit.c

LIBOBJ = $(CXXSRC:%.cc=%.o) $(CSRC:%.c=%.o) $(GENCXXSRC:%.cc=%.o)

INC = volume.h shape.h volume_i.h implicit.h

DEPENDINCLUDE = $(GENINC) $(INC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = isosurf

DISTFILES = isotest.cc $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h

default:: $(DEPENDINCLUDE)
default:: isosurf.a

isotest: isotest.o isosurf.a util/options/liboptions.a \
	$(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LD) $(LDFLAGS) -o isotest $^ $(SYSLIBS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d) isotest.d: $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) isotest.d
endif

####### Rules for generating source code ######################################

edgeMap.cc edgeMap.h:
	$(GENCLASS) -2 RefVertex ref RefEdgeAVLSet ref Map edge
	mv edgeMap.h tmp.edgeMap.h
	grep -v edgedefs.h tmp.edgeMap.h > edgeMap.h
	rm tmp.edgeMap.h
	prepend-header edgeMap.h "#include <math/isosurf/edgeAVLSet.h>"

edgeRAVLMap.cc edgeRAVLMap.h:
	$(GENCLASS) -2 RefVertex ref RefEdgeAVLSet ref RAVLMap edge

edgeSet.cc edgeSet.h:
	$(GENCLASS) RefEdge ref Set edge
	mv edgeSet.h tmp.edgeSet.h
	grep -v edgedefs.h tmp.edgeSet.h > edgeSet.h
	rm tmp.edgeSet.h
	prepend-header edgeSet.h "#include <math/isosurf/edge.h>"
	prepend-header edgeSet.h "#include <math/isosurf/edgedefs.h>"

edgeAVLSet.cc edgeAVLSet.h:
	$(GENCLASS) RefEdge ref AVLSet edge
