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
VPATH+=$(SRCDIR)/levmar

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

LD=$(CXX)

CXXSRC = function.cc opt.cc update.cc qnewton.cc newton.cc \
         gdiis.cc efc.cc steep.cc dfp.cc powell.cc \
         scextrap.cc scextrapmat.cc diis.cc transform.cc conv.cc
CSRC = lm.c Axb.c misc.c lmlec.c lmbc.c

LIBOBJ = $(CXXSRC:%.cc=%.$(OBJSUF)) $(CSRC:%.c=%.$(OBJSUF))

INC = function.h opt.h gdiis.h efc.h update.h qnewton.h newton.h \
      scextrap.h scextrapmat.h diis.h transform.h conv.h \
      lm.h

DEPENDINCLUDE =

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libSCoptimize

kvopt.testrun: DO_TESTRUN=no

TESTPROGS = kvopt opttest scextest lmdemo

DISTFILES = opttest.cc $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

default:: $(DEPENDINCLUDE)

install_inc::
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/levmar/lm.h $(installroot)$(includedir)/$(SRCDIR:$(SRCTOPDIR)/src/lib/%=%)/levmar

opttest: opttest.$(OBJSUF) $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LTLINK) $(LD) $(LDFLAGS) -o opttest $^ $(SYSLIBS) $(LTLINKBINOPTS)

opttest.$(OBJSUF): opttest.cc
	$(LTCOMP) $(CXX) $(CXXFLAGS) $(CPPFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

scextest: scextest.$(OBJSUF) libSCoptimize.$(LIBSUF) \
	$(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LTLINK) $(LD) $(LDFLAGS) -o scextest $^ $(SYSLIBS) $(LTLINKBINOPTS)

scextest.$(OBJSUF): scextest.cc
	$(LTCOMP) $(CXX) $(CXXFLAGS) $(CPPFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

kvopt: kvopt.$(OBJSUF) $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LTLINK) $(LD) $(LDFLAGS) -o kvopt $^ $(SYSLIBS) $(LTLINKBINOPTS)

lmdemo: lmdemo.$(OBJSUF) $(LIBS)
	$(LTLINK) $(LD) $(LDFLAGS) -o lmdemo $^ $(SYSLIBS) $(LTLINKBINOPTS)

Axb.$(OBJSUF) Axb.d: f77sym.h
lmlec.$(OBJSUF) lmlec.d: f77sym.h
misc.$(OBJSUF) misc.d: f77sym.h

f77sym.h: f77sym.in
	$(MKF77SYM) $< $@

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

kvopt.d opttest.d scextest.d lmdemo.d $(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) kvopt.d opttest.d scextest.d lmdemo.d
endif

