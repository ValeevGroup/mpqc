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

TARGET_TO_MAKE = libmisc
BIN_OR_LIB = LIB

CSRC = ieee.c checkalloc.c mtype.c

CXXSRC = compute.cc formio.cc bug.cc regtime.cc exenv.cc

TESTSRC = comptest.cc formiot.cc bugtest.cc tregtime.cc

INC = libmisc.h memory.h timer.h compute.h ieee.h bug.h exenv.h

DISTFILES = $(CSRC) $(CXXSRC) $(TESTSRC) $(INC) Makefile LIBS.h

LIBOBJ= $(CSRC:%.c=%.$(OBJSUF)) \
        $(CXXSRC:%.cc=%.$(OBJSUF))

TESTOBJ= $(TESTSRC:%.cc=%.$(OBJSUF))

DEPENDINCLUDE = $(INC)

#####################################################################

default:: $(DEPENDINCLUDE)

comptest: comptest.cc $(TOPDIR)/lib/libmisc.a $(TOPDIR)/lib/libstate.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o comptest

tregtime: tregtime.cc $(TOPDIR)/lib/libmisc.a $(TOPDIR)/lib/libref.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o tregtime

bugtest: bugtest.cc \
	$(TOPDIR)/lib/libkeyval.a $(TOPDIR)/lib/libstate.a \
	$(TOPDIR)/lib/libclass.a $(TOPDIR)/lib/libmisc.a \
        $(TOPDIR)/lib/libgroup.a \
	$(TOPDIR)/lib/libcontainer.a $(TOPDIR)/lib/libref.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" $^ $(SYSLIBS) -o bugtest

formiot: formiot.cc libmisc.a libgroup.a libkeyval.a libstate.a libclass.a \
        libcontainer.a libref.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o formiot $(SYSLIBS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d) $(TESTOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTOBJ:.o=.d)
endif

