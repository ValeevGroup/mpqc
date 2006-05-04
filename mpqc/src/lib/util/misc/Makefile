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
include $(TOPDIR)/lib/Makedirlist

# needed for corba
CXXFLAGS := $(CXXFLAGS:-fno-implicit-templates=)

TARGET_TO_MAKE = libSCmisc
BIN_OR_LIB = LIB

CSRC =

CXXSRC = compute.cc formio.cc bug.cc regtime.cc exenv.cc units.cc ieee.cc

TESTSRC = comptest.cc formiot.cc bugtest.cc tregtime.cc unittest.cc \
          autovectest.cc

TESTPROGS = autovectest comptest formiot bugtest tregtime unittest

INC = timer.h compute.h ieee.h bug.h exenv.h

DISTFILES = $(CSRC) $(CXXSRC) $(TESTSRC) $(INC) Makefile LIBS.h

LIBOBJ= $(CSRC:%.c=%.$(OBJSUF)) \
        $(CXXSRC:%.cc=%.$(OBJSUF))

TESTOBJ= $(TESTSRC:%.cc=%.$(OBJSUF))

DEPENDINCLUDE = $(INC)

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_CCA),yes)
  CXXSRC += ccaenv.cc
  INC += ccaenv.h
  LTLINKLIBOPTS += -L$(CCAFE_LIB) -R$(CCAFE_LIB) -lccafeCore \
    -L$(BABEL_LIB) -R$(BABEL_LIB) -lsidl
  CPPFLAGS += -I../../chemistry/cca
endif

#####################################################################

default:: $(DEPENDINCLUDE)

autovectest: autovectest.$(OBJSUF)
	$(LTLINK) $(LD) $(LDFLAGS) $^ -o $@ $(LTLINKBINOPTS)

comptest: comptest.$(OBJSUF) $(LIBS)
	$(LTLINK) $(LD) $(LDFLAGS) $^ -o comptest $(SYSLIBS) $(LTLINKBINOPTS)

tregtime: tregtime.$(OBJSUF) $(LIBS)
	$(LTLINK) $(LD) $(LDFLAGS) $^ -o tregtime $(SYSLIBS) $(LTLINKBINOPTS)

scinttest: scinttest.$(OBJSUF) $(LIBS)
	$(LTLINK) $(LD) $(LDFLAGS) $^ -o scinttest $(SYSLIBS) $(LTLINKBINOPTS)

bugtest: bugtest.$(OBJSUF) $(LIBS)
	$(LTLINK) $(LD) $(LDFLAGS) $^ $(SYSLIBS) -o bugtest $(LTLINKBINOPTS)

scextest: scextest.$(OBJSUF) $(LIBS)
	$(LTLINK) $(LD) $(LDFLAGS) $^ $(SYSLIBS) -o $@ $(LTLINKBINOPTS)

bugtest.$(OBJSUF): bugtest.cc
	$(LTCOMP) $(CXX) $(CXXFLAGS) $(CPPFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $< $(LTLINKBINOPTS)

unittest: unittest.$(OBJSUF) $(LIBS) 
	$(LTLINK) $(LD) $(LDFLAGS) $^ $(SYSLIBS) -o unittest $(LTLINKBINOPTS)

formiot: formiot.$(OBJSUF) $(LIBS)
	$(LTLINK) $(LD) $(LDFLAGS) $^ -o formiot $(SYSLIBS) $(LTLINKBINOPTS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

distclean::
	/bin/rm -f autovectest comptest tregtime bugtest formiot unittest

$(LIBOBJ:.$(OBJSUF)=.d) $(TESTOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) $(TESTOBJ:.$(OBJSUF)=.d)
endif

# need babel headers generated
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_CCA),yes)
  ifneq ($(DODEPEND),no)
    include ../../chemistry/cca/stamp-babel
  endif
endif
../../chemistry/cca/stamp-babel:
	cd ../../chemistry/cca; $(MAKE) stamp-babel
../../../../lib/cca/libSCchemistry.$(LIBSUF):
	cd ../../chemistry/cca; $(MAKE)

