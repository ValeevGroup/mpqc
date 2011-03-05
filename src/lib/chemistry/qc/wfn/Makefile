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

TARGET_TO_MAKE = libSCwfn
BIN_OR_LIB = LIB

TESTPROGS = wfntest

SRCS = wfn.cc densval.cc natorbval.cc obwfn.cc density.cc hcorewfn.cc \
       accum.cc nao.cc orbital.cc solvent.cc eht.cc esp.cc femo.cc \
       orbitalspace.cc orbitalspace_utils.cc rdm.cc spin.cc

INC = wfn.h accum.h density.h eht.h esp.h femo.h obwfn.h orbital.h \
      solvent.h \
      orbitalspace.h orbitalspace.timpl.h \
      orbitalspace_utils.h orbitalspace_utils.timpl.h rdm.h \
      spin.h

LIBOBJ= $(SRCS:%.cc=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

WFNTESTLIBS: = $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/../dft/LIBS.h)

wfntest:: wfntest.$(OBJSUF) $(WFNTESTLIBS)
	$(LTLINK) $(CXX) $(LDFLAGS) -o wfntest $^ $(SYSLIBS) $(LTLINKBINOPTS)

wfntest.$(OBJSUF): wfntest.cc
	$(LTCOMP) $(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

$(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) $(TESTPROGS:%=%.d)
endif

