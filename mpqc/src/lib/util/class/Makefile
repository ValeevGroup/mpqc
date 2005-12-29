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

CXXSRC = class.cc scexception.cc

LIBOBJ = $(CXXSRC:%.cc=%.$(OBJSUF)) $(GENCXXSRC:%.cc=%.$(OBJSUF))

INC = class.h proxy.h scexception.h

DEPENDINCLUDE = $(INC) $(GENINC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libSCclass

TESTCXXSRC = classtest.cc scextest.cc
DISTFILES = $(CXXSRC) $(INC) Makefile LIBS.h $(TESTCXXSRC)

TESTPROGS = classtest scextest

default:: $(DEPENDINCLUDE)

interface:: $(DEPENDINCLUDE)

# a test program
LD = $(CXX)
classtest: classtest.$(OBJSUF) libSCclass.$(LIBSUF) libSCcontainer.$(LIBSUF) libSCref.$(LIBSUF) libSCmisc.$(LIBSUF) libSCstate.$(LIBSUF) libSCkeyval.$(LIBSUF)
	$(LTLINK) $(LD) $(LDFLAGS) -o classtest $^ $(SYSLIBS) $(LTLINKBINOPTS)

scextest: scextest.$(OBJSUF) libSCclass.$(LIBSUF) libSCcontainer.$(LIBSUF) libSCref.$(LIBSUF) libSCmisc.$(LIBSUF) libSCstate.$(LIBSUF) libSCkeyval.$(LIBSUF)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

bug: bug.cc
	$(LTLINK) $(CXX) $(CXXFLAGS) -o bug $^ $(LTLINKBINOPTS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(TESTCXXSRC:%.cc=%.d) $(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(TESTCXXSRC:%.cc=%.d) $(LIBOBJ:.$(OBJSUF)=.d)
endif
