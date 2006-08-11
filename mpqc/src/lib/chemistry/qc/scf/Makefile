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
TOPDIR=../../../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

include $(TOPDIR)/lib/Makedirlist

TARGET_TO_MAKE = libSCscf
BIN_OR_LIB = LIB

TESTPROGS = scftest focktest

SRCS = fbclhf.cc scf.cc scfvector.cc scfgradient.cc scfops.cc effh.cc \
       clscf.cc clhf.cc hsosscf.cc hsoshf.cc \
       ossscf.cc osshf.cc tcscf.cc tchf.cc \
       uscf.cc uhf.cc fockbuild.cc clhfcontrib.cc \
       tmplinst.cc svd.cc

LIBOBJ= $(SRCS:%.cc=%.$(OBJSUF))

#Don't generate library dependencies for scf
#LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

SCFTESTLIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/../dft/LIBS.h)

scftest:: scftest.$(OBJSUF) $(SCFTESTLIBS)
	$(LTLINK) $(CXX) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

scftest.$(OBJSUF): scftest.cc
	$(LTCOMP) $(CXX) $(CXXFLAGS) $(CPPFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

focktest:: focktest.$(OBJSUF) $(SCFTESTLIBS)
	$(LTLINK) $(CXX) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

focktest.$(OBJSUF): focktest.cc
	$(LTCOMP) $(CXX) $(CXXFLAGS) $(CPPFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

$(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) scftest.d focktest.d
endif

