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

TOPDIR=../../../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libSCintv3
BIN_OR_LIB = LIB

TESTCSRC = inttest.cc fjttest.cc

CXXSRC = offsets.cc tformv3.cc \
         int1e.cc comp1e.cc \
         int2e.cc comp2e.cc build2e.cc print2e.cc shift2e.cc init2e.cc \
         comp2e3c.cc \
         bounds.cc store.cc storage.cc \
         obintv3.cc tbintv3.cc intv3.cc \
         array.cc

INC = int1e.h int2e.h utils.h tformv3.h \
      tbintv3.h obintv3.h intv3.h cartitv3.h

LIBSRC = $(CXXSRC) $(CSRC) $(XCSRC)
LIBOBJ = $(LIBSRC:%.c=%.$(OBJSUF))
LIBOBJ := $(LIBOBJ:%.cc=%.$(OBJSUF))

DISTFILES = $(TESTCSRC) $(INC) atoms.sgen Makefile LIBS.h

DEPENDINCLUDE = $(INC) $(GENINC) $(SGENINC)

TESTPROGS = inttest fjttest

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

#############################################################

default:: $(DEPENDINCLUDE)
interface:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

fjttest: fjttest.$(OBJSUF) \
	$(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LTLINK) $(CXX) $(LDFLAGS) -o fjttest $^ $(SYSLIBS) $(LTLINKBINOPTS)

inttest: inttest.$(OBJSUF) \
	$(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LTLINK) $(CXX) $(LDFLAGS) -o inttest $^ $(SYSLIBS) $(LTLINKBINOPTS)

inttest.$(OBJSUF): inttest.cc
	$(LTCOMP) $(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

$(TESTOBJ:.$(OBJSUF)=.d) $(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) $(TESTOBJ:.$(OBJSUF)=.d)
endif

#############################################################
