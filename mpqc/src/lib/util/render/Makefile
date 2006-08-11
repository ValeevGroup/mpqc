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
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

CXXSRC = render.cc oogl.cc object.cc polygons.cc polysphere.cc \
         transform.cc appearance.cc material.cc sphere.cc algebra3.cc \
         color.cc polylines.cc tempinst.cc animate.cc

LIBOBJ = $(CXXSRC:%.cc=%.$(OBJSUF))

INC = render.h oogl.h object.h sphere.h color.h \
      algebra3.h parameter.h find.h stack.h polygons.h polysphere.h \
      polylines.h animate.h

DEPENDINCLUDE = $(INC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libSCrender

TESTSRC = rentest.cc
TESTOBJ = $(TESTSRC:%.cc=%.$(OBJSUF))
TESTFILES = 

TESTPROGS = rentest

DISTFILES = $(CXXSRC) $(INC) Makefile $(TESTSRC) $(TESTFILES)

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

default:: $(DEPENDINCLUDE)

interface:: $(DEPENDINCLUDE)

LD = $(CXX)
rentest: $(TESTOBJ) \
	 $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LTLINK) $(LD) $(LDFLAGS) -o rentest $^ $(SYSLIBS) $(LTLINKBINOPTS)

rentest.$(OBJSUF): rentest.cc
	$(LTCOMP) $(CXX) -DSRCDIR=\"$(SRCDIR)\" $(CPPFLAGS) $(CXXFLAGS) -c $<

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) $(TESTOBJ:%.$(OBJSUF)=%.d)
endif
