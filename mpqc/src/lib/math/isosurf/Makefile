#
# Makefile
#
# Copyright (C) 1996 Limit Point Systems, Inc.
#
# Author: Curtis Janssen <cljanss@ca.sandia.gov>
# Maintainer: LPS
#
# This file is part of the SC Toolkit.
#
# The SC Toolkit is free software; you can redistribute it and/or modify
# it under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# The SC Toolkit is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General Public License
# along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
# the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
#
# The U.S. Government is granted a limited license as per AL 91-7.
#

TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

LD=$(CXX)

GENSRC = edgeSet.cc edgeAVLSet.cc \
            edgeMap.cc edgeRAVLMap.cc \
            vertexSet.cc vertexAVLSet.cc \
            veMap.cc veRAVLMap.cc \
            vtsMap.cc vtsRAVLMap.cc \
            triSet.cc triAVLSet.cc \
            triMap.cc triRAVLMap.cc
GENINC = $(GENSRC:%.cc=%.h)
CXXSRC = volume.cc shape.cc volume_i.cc \
         surf.cc surfse.cc surfor.cc surfst.cc \
         vertex.cc edge.cc triangle.cc isosurf.cc \
         tricoef.cc
CSRC = implicit.c

CLASSES = Volume Shape UnionShape SphereShape \
          UncappedTorusHoleShape NonreentrantUncappedTorusHoleShape \
          ReentrantUncappedTorusHoleShape Uncapped5SphereExclusionShape

LIBOBJ = $(CXXSRC:%.cc=%.o) $(CSRC:%.c=%.o) $(GENSRC:%.cc=%.o)

INC = volume.h shape.h volume_i.h implicit.h tricoef.h

DEPENDINCLUDE = $(GENINC) $(INC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libisosurf

DISTFILES = isotest.cc $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h

default:: $(DEPENDINCLUDE)

isotest.o: isotest.cc
	$(CXX) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

isotest: isotest.o \
	$(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LD) $(LDFLAGS) -o isotest $^ $(SYSLIBS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d) isotest.d: $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) isotest.d
endif

####### Rules for generating source code ######################################

#tricoefMap.cc tricoefMap.h:
#	$(GENCLASS) -2 TriInterpCoefKey ref RefTriInterpCoef ref Map tricoef
#	mv tricoefMap.h tmp.tricoefMap.h
#	grep -v tricoefdefs.h tmp.tricoefMap.h > tricoefMap.h
#	rm tmp.tricoefMap.h
#	$(PREPEND_HEADER) tricoefMap.h "#include <math/isosurf/tricoef.h>"
#
#tricoefRAVLMap.cc tricoefRAVLMap.h:
#	$(GENCLASS) -2 TriInterpCoefKey ref RefTriInterpCoef ref RAVLMap tricoef

#######

veMap.cc veMap.h:
	$(GENCLASS) -2 RefVertex ref RefEdge ref Map ve
	mv veMap.h tmp.veMap.h
	grep -v vedefs.h tmp.veMap.h > veMap.h
	rm tmp.veMap.h
	$(PREPEND_HEADER) veMap.h "#include <math/isosurf/vertex.h>"
	$(PREPEND_HEADER) veMap.h "#include <math/isosurf/edge.h>"
	$(PREPEND_HEADER) veMap.h "#include <math/isosurf/edgedefs.h>"

veRAVLMap.cc veRAVLMap.h:
	$(GENCLASS) -2 RefVertex ref RefEdge ref RAVLMap ve

#######

vtMap.cc vtMap.h:
	$(GENCLASS) -2 RefVertex ref RefTriangle ref Map vt
	mv vtMap.h tmp.vtMap.h
	grep -v vtdefs.h tmp.vtMap.h > vtMap.h
	rm tmp.vtMap.h
	$(PREPEND_HEADER) vtMap.h "#include <math/isosurf/vertex.h>"
	$(PREPEND_HEADER) vtMap.h "#include <math/isosurf/triangle.h>"
	$(PREPEND_HEADER) vtMap.h "#include <math/isosurf/edgedefs.h>"

vtRAVLMap.cc vtRAVLMap.h:
	$(GENCLASS) -2 RefVertex ref RefTriangle ref RAVLMap vt

#######

vtsMap.cc vtsMap.h:
	$(GENCLASS) -2 RefVertex ref RefTriangleAVLSet ref Map vts
	mv vtsMap.h tmp.vtsMap.h
	grep -v vtsdefs.h tmp.vtsMap.h > vtsMap.h
	rm tmp.vtsMap.h
	$(PREPEND_HEADER) vtsMap.h "#include <math/isosurf/vertex.h>"
	$(PREPEND_HEADER) vtsMap.h "#include <math/isosurf/triAVLSet.h>"
	$(PREPEND_HEADER) vtsMap.h "#include <math/isosurf/edgedefs.h>"

vtsRAVLMap.cc vtsRAVLMap.h:
	$(GENCLASS) -2 RefVertex ref RefTriangleAVLSet ref RAVLMap vts

#######

vertexSet.cc vertexSet.h:
	$(GENCLASS) RefVertex ref Set vertex
	mv vertexSet.h tmp.vertexSet.h
	grep -v vertexdefs.h tmp.vertexSet.h > vertexSet.h
	rm tmp.vertexSet.h
	$(PREPEND_HEADER) vertexSet.h "#include <math/isosurf/vertex.h>"
	$(PREPEND_HEADER) vertexSet.h "#include <math/isosurf/edgedefs.h>"

vertexAVLSet.cc vertexAVLSet.h:
	$(GENCLASS) RefVertex ref AVLSet vertex

#######

edgeMap.cc edgeMap.h:
	$(GENCLASS) -2 RefVertex ref RefEdgeAVLSet ref Map edge
	mv edgeMap.h tmp.edgeMap.h
	grep -v edgedefs.h tmp.edgeMap.h > edgeMap.h
	rm tmp.edgeMap.h
	$(PREPEND_HEADER) edgeMap.h "#include <math/isosurf/edgeAVLSet.h>"

edgeRAVLMap.cc edgeRAVLMap.h:
	$(GENCLASS) -2 RefVertex ref RefEdgeAVLSet ref RAVLMap edge

edgeSet.cc edgeSet.h:
	$(GENCLASS) RefEdge ref Set edge
	mv edgeSet.h tmp.edgeSet.h
	grep -v edgedefs.h tmp.edgeSet.h > edgeSet.h
	rm tmp.edgeSet.h
	$(PREPEND_HEADER) edgeSet.h "#include <math/isosurf/edge.h>"
	$(PREPEND_HEADER) edgeSet.h "#include <math/isosurf/edgedefs.h>"

edgeAVLSet.cc edgeAVLSet.h:
	$(GENCLASS) RefEdge ref AVLSet edge

###

triMap.cc triMap.h:
	$(GENCLASS) -2 RefEdge ref RefTriangleAVLSet ref Map tri
	mv triMap.h tmp.triMap.h
	grep -v tridefs.h tmp.triMap.h > triMap.h
	rm tmp.triMap.h
	$(PREPEND_HEADER) triMap.h "#include <math/isosurf/triAVLSet.h>"

triRAVLMap.cc triRAVLMap.h:
	$(GENCLASS) -2 RefEdge ref RefTriangleAVLSet ref RAVLMap tri

triSet.cc triSet.h:
	$(GENCLASS) RefTriangle ref Set tri
	mv triSet.h tmp.triSet.h
	grep -v tridefs.h tmp.triSet.h > triSet.h
	rm tmp.triSet.h
	$(PREPEND_HEADER) triSet.h "#include <math/isosurf/triangle.h>"
	$(PREPEND_HEADER) triSet.h "#include <math/isosurf/edgedefs.h>"

triAVLSet.cc triAVLSet.h:
	$(GENCLASS) RefTriangle ref AVLSet tri
