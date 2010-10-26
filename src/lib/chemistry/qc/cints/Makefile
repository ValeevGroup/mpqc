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

TARGET_TO_MAKE = libSCcints
BIN_OR_LIB = LIB

#DEFINES += -DDMALLOC

TESTCSRC = cintstest.cc

CXXSRC = cints.cc obintcints.cc tbintcints.cc tform.cc int2e.cc eri.cc comp_eri.cc \
        storage.cc grt.cc comp_grt.cc \
	primpairs.cc shellpairs.cc permute2e.cc int1e.cc overlap.cc kinetic.cc nuclear.cc \
	hcore.cc obosrr.cc edipole.cc equadrupole.cc

INC = cints.h obintcints.h tbintcints.h tform.h cartit.h macros.h primpairs.h shellpairs.h

LIBSRC = $(CXXSRC) $(CSRC) $(XCSRC)
LIBOBJ = $(LIBSRC:%.c=%.$(OBJSUF))
LIBOBJ := $(LIBOBJ:%.cc=%.$(OBJSUF))

LTLINKLIBOPTS += $(LDFLAGS)
ifeq ($(HAVE_LIBINT),yes)
  LTLINKLIBOPTS += -lint
endif
ifeq ($(HAVE_LIBR12),yes)
  LTLINKLIBOPTS += -lr12
endif
ifeq ($(HAVE_LIBDERIV),yes)
  LTLINKLIBOPTS += -lderiv
endif

DISTFILES = $(TESTCSRC) $(INC) Makefile

DEPENDINCLUDE = $(INC) $(GENINC) $(SGENINC)

#############################################################

default:: $(DEPENDINCLUDE)
interface:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

cintstest: cintstest.$(OBJSUF) \
	$(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LTLINK) $(CXX) $(LDFLAGS) -o cintstest $^ $(SYSLIBS) $(LTLINKBINOPTS)

cintstest.$(OBJSUF): cintstest.cc
	$(LTCOMP) $(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

targetclean::
	-rm -f cintstest

$(TESTOBJ:.$(OBJSUF)=.d) $(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) $(TESTOBJ:.$(OBJSUF)=.d)
endif

#############################################################
