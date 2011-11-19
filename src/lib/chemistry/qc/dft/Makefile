#
# Makefile
#
# Copyright (C) 1997 Limit Point Systems, Inc.
#
# Author: Curtis Janssen <cljanss@limitpt.com>
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

TARGET_TO_MAKE = libSCdft
BIN_OR_LIB = LIB

TESTPROGS = dfttest lebedev

SRCS = functional.cc integrator.cc \
       clks.cc uks.cc hsosks.cc lebedev.c \
       tmplinst.cc am05.cc solvent.cc

INC = functional.h integrator.h \
      clks.h uks.h hsosks.h lebedev.h \
      am05.h solvent.h

LIBOBJ= $(SRCS:%.cc=%.$(OBJSUF))
LIBOBJ:= $(LIBOBJ:%.c=%.$(OBJSUF))
LIBOBJ:= $(LIBOBJ:%.f=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

dfttest:: dfttest.$(OBJSUF) $(LIBS)
	$(LTLINK) $(CXX) $(LDFLAGS) -o dfttest $^ $(SYSLIBS) $(LTLINKBINOPTS)

dfttest.$(OBJSUF): dfttest.cc
	$(LTCOMP) $(CXX) $(CXXFLAGS) $(CPPFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

lebedev: lebedev.c
	$(LTLINK) $(CC) $(CFLAGS) $(CPPFLAGS) -DDEFINE_MAIN -o $@ $< $(LTLINKBINOPTS) -lm

DEPFILES:=$(LIBOBJ:.$(OBJSUF)=.d)

$(DEPFILES): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(DEPFILES) dfttest.d
endif
