#
# Makefile
#
# Copyright (C) 1996 Limit Point Systems, Inc.
#
# Author: Edward Seidl <seidl@janed.com>
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

DEFINES += -DSRCDIR=\"$(SRCDIR)\"

TESTPROGS = test_qcfnc1

TARGET_TO_MAKE = libSCexternqc
BIN_OR_LIB = LIB

CXXSRC = qc_fnc.cc

INC = qc_fnc.h

DISTFILES = $(CXXSRC) $(INC) Makefile $(TESTSRC) $(TESTFILES)

LIBOBJ= $(CXXSRC:%.cc=%.$(OBJSUF))

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) $(TESTOBJ:%.$(OBJSUF)=%.d)
endif

TESTSRC += test_qcfnc0.cc
test_qcfnc0: test_qcfnc0.$(OBJSUF) $(LIBS) 
	$(LTLINK) $(CXX) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

test_qcfnc0.$(OBJSUF): test_qcfnc0.cc
	$(LTCOMP) $(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $<

TESTSRC += test_qcfnc1.c
test_qcfnc1: test_qcfnc1.$(OBJSUF) $(LIBS) 
	$(LTLINK) $(CXX) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

test_qcfnc1.$(OBJSUF): test_qcfnc1.c
	$(LTCOMP) $(CC) $(CCFLAGS) $(CPPFLAGS) -c $<

TESTSRC += test_qcfnc2.f90
test_qcfnc2: test_qcfnc2.$(OBJSUF) $(LIBS) 
	$(LTLINK) $(CXX) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

test_qcfnc2.$(OBJSUF): test_qcfnc2.f90
	$(LTCOMP) $(FC) $(FFLAGS) $(CPPFLAGS) -c $<
	