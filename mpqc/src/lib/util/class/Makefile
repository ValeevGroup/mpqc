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

CXXSRC = class.cc proxy.cc

LIBOBJ = $(CXXSRC:%.cc=%.o) $(GENCXXSRC:%.cc=%.o)

CLASSES = DescribedClass DescribedClassProxy

INC = class.h proxy.h
GENINC = clasmacr.h

DEPENDINCLUDE = $(INC) $(GENINC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libclass

TESTCXXSRC = classtest.cc
DISTFILES = $(CXXSRC) $(INC) Makefile LIBS.h $(TESTCXXSRC)

default:: $(DEPENDINCLUDE)
default:: dlinfo

interface:: $(DEPENDINCLUDE)

clasmacr.h: clastmpl.h $(TMPL2MACDEP)
	$(TMPL2MAC) $< > clasmacr.h

# a test program
LD = $(CXX)
classtest: classtest.o libclass.$(LIBSUF) libcontainer.$(LIBSUF) libref.$(LIBSUF) libmisc.$(LIBSUF) libstate.$(LIBSUF) libkeyval.$(LIBSUF)
	$(LD) $(LDFLAGS) -o classtest $^ $(SYSLIBS)

bug: bug.cc
	$(CXX) $(CXXFLAGS) -o bug $^

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(TESTCXXSRC:%.cc=%.d) $(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(TESTCXXSRC:%.cc=%.d) $(LIBOBJ:.o=.d)
endif
