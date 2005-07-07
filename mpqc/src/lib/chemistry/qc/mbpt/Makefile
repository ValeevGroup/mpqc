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

TARGET_TO_MAKE = libSCmbpt
BIN_OR_LIB = LIB

TESTPROGS = mbpttest

CXXSRCS = mbpt.cc csgrad.cc csgrade12.cc csgrads2pdm.cc csgmat.cc cscphf.cc \
          hsosv1.cc hsosv1e1.cc hsosv2.cc hsosv2lb.cc util.cc mp2extrap.cc \
          csgrad34qb.cc

CSRCS = bzerofast.c

LIBOBJ= $(CXXSRCS:%.cc=%.$(OBJSUF)) $(CSRCS:%.c=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d)
endif

LIBS = $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

mbpttest:: mbpttest.$(OBJSUF) $(LIBS) libSCdft.$(LIBSUF)
	$(LTLINK) $(CXX) $(LDFLAGS) -o mbpttest $^ $(SYSLIBS) $(LTLINKBINOPTS)

mbpttest.$(OBJSUF): mbpttest.cc
	$(LTCOMP) $(CXX) $(CXXFLAGS) $(CPPFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<
