TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

LD=$(CXX)

GENCXXSRC = edgeSet.cc edgeAVLSet.cc \
            edgeMap.cc edgeRAVLMap.cc \
            vertexSet.cc vertexAVLSet.cc \
            veMap.cc veRAVLMap.cc \
            vtsMap.cc vtsRAVLMap.cc \
            triSet.cc triAVLSet.cc \
            triMap.cc triRAVLMap.cc
GENINC = $(GENCXXSRC:%.cc=%.h)
CXXSRC = volume.cc shape.cc volume_i.cc \
         surf.cc surfse.cc surfor.cc surfst.cc \
         vertex.cc edge.cc triangle.cc isosurf.cc
CSRC = implicit.c

LIBOBJ = $(CXXSRC:%.cc=%.o) $(CSRC:%.c=%.o) $(GENCXXSRC:%.cc=%.o)

INC = volume.h shape.h volume_i.h implicit.h

DEPENDINCLUDE = $(GENINC) $(INC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = isosurf

DISTFILES = isotest.cc $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h

default:: $(DEPENDINCLUDE)
default:: $(TOPDIR)/lib/isosurf.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/isosurf.so
endif

isotest: isotest.o isosurf.a util/options/liboptions.a \
	$(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LD) $(LDFLAGS) -o isotest $^ $(SYSLIBS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d) isotest.d: $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) isotest.d
endif

####### Rules for generating source code ######################################

veMap.cc veMap.h:
	$(GENCLASS) -2 RefVertex ref RefEdge ref Map ve
	mv veMap.h tmp.veMap.h
	grep -v vedefs.h tmp.veMap.h > veMap.h
	rm tmp.veMap.h
	prepend-header veMap.h "#include <math/isosurf/vertex.h>"
	prepend-header veMap.h "#include <math/isosurf/edge.h>"
	prepend-header veMap.h "#include <math/isosurf/edgedefs.h>"

veRAVLMap.cc veRAVLMap.h:
	$(GENCLASS) -2 RefVertex ref RefEdge ref RAVLMap ve

#######

vtMap.cc vtMap.h:
	$(GENCLASS) -2 RefVertex ref RefTriangle ref Map vt
	mv vtMap.h tmp.vtMap.h
	grep -v vtdefs.h tmp.vtMap.h > vtMap.h
	rm tmp.vtMap.h
	prepend-header vtMap.h "#include <math/isosurf/vertex.h>"
	prepend-header vtMap.h "#include <math/isosurf/triangle.h>"
	prepend-header vtMap.h "#include <math/isosurf/edgedefs.h>"

vtRAVLMap.cc vtRAVLMap.h:
	$(GENCLASS) -2 RefVertex ref RefTriangle ref RAVLMap vt

#######

vtsMap.cc vtsMap.h:
	$(GENCLASS) -2 RefVertex ref RefTriangleAVLSet ref Map vts
	mv vtsMap.h tmp.vtsMap.h
	grep -v vtsdefs.h tmp.vtsMap.h > vtsMap.h
	rm tmp.vtsMap.h
	prepend-header vtsMap.h "#include <math/isosurf/vertex.h>"
	prepend-header vtsMap.h "#include <math/isosurf/triAVLSet.h>"
	prepend-header vtsMap.h "#include <math/isosurf/edgedefs.h>"

vtsRAVLMap.cc vtsRAVLMap.h:
	$(GENCLASS) -2 RefVertex ref RefTriangleAVLSet ref RAVLMap vts

#######

vertexSet.cc vertexSet.h:
	$(GENCLASS) RefVertex ref Set vertex
	mv vertexSet.h tmp.vertexSet.h
	grep -v vertexdefs.h tmp.vertexSet.h > vertexSet.h
	rm tmp.vertexSet.h
	prepend-header vertexSet.h "#include <math/isosurf/vertex.h>"
	prepend-header vertexSet.h "#include <math/isosurf/edgedefs.h>"

vertexAVLSet.cc vertexAVLSet.h:
	$(GENCLASS) RefVertex ref AVLSet vertex

#######

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

###

triMap.cc triMap.h:
	$(GENCLASS) -2 RefEdge ref RefTriangleAVLSet ref Map tri
	mv triMap.h tmp.triMap.h
	grep -v tridefs.h tmp.triMap.h > triMap.h
	rm tmp.triMap.h
	prepend-header triMap.h "#include <math/isosurf/triAVLSet.h>"

triRAVLMap.cc triRAVLMap.h:
	$(GENCLASS) -2 RefEdge ref RefTriangleAVLSet ref RAVLMap tri

triSet.cc triSet.h:
	$(GENCLASS) RefTriangle ref Set tri
	mv triSet.h tmp.triSet.h
	grep -v tridefs.h tmp.triSet.h > triSet.h
	rm tmp.triSet.h
	prepend-header triSet.h "#include <math/isosurf/triangle.h>"
	prepend-header triSet.h "#include <math/isosurf/edgedefs.h>"

triAVLSet.cc triAVLSet.h:
	$(GENCLASS) RefTriangle ref AVLSet tri
