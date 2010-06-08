#
# Makefile
#
# Copyright (C) 1996 Limit Point Systems, Inc.
#
# Author: Ida Nielsen <ibniels@ca.sandia.gov>
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

TARGET_TO_MAKE = libSClmp2
BIN_OR_LIB = LIB

TESTPROGS = testsma timesma avlmmaptest

CXXSRCS = collective.cc dgemminfo.cc domain.cc extrap.cc lcorr.cc lmp2.cc \
          pack2e.cc pop_local.cc sma.cc smautil.cc

LIBOBJ= $(CXXSRCS:%.cc=%.$(OBJSUF)) $(CSRCS:%.c=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d)
endif

domain.$(OBJSUF) domain.d smautil.d : f77sym.h

f77sym.h: f77sym.h.in
	$(MKF77SYM) $< $@

distclean::
	/bin/rm -f f77sym.h

TESTSRC += testsma.cc
LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)
testsma: testsma.o $(SMAOBJ)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $^ -L$(LIBDIR) $(LIBS)

TESTSRC += timesma.cc
timesma: timesma.o $(SMAOBJ)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $^ -L$(LIBDIR) $(LIBS)

TESTSRC += avlmmaptest.cc
avlmmaptest: avlmmaptest.o
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $^ -L$(LIBDIR) $(LIBS)
