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

TARGET_TO_MAKE = libSClibint2
BIN_OR_LIB = LIB

#DEFINES += -DDMALLOC

TESTCSRC = libint2test.cc

CXXSRC = libint2.cc obintlibint2.cc tbintlibint2.cc tform.cc int2e.cc bounds.cc storage.cc \
	primpairs.cc shellpairs.cc permute2e.cc int1e.cc overlap.cc kinetic.cc nuclear.cc \
	hcore.cc obosrr.cc edipole.cc equadrupole.cc p4.cc gto.cc
ifeq ($(HAVE_LIBINT2_ERI),yes)
  CXXSRC += eri.cc comp_eri.cc
  CXXSRC += g12nc.cc comp_g12nc.cc
endif
ifeq ($(HAVE_LIBINT2_G12),yes)
  ifeq ($(HAVE_LIBINT2_G12_T1G12),yes)
    CXXSRC += g12.cc comp_g12.cc
  endif
endif
ifeq ($(HAVE_LIBINT2_G12DKH),yes)
  CXXSRC += g12dkh.cc comp_g12dkh.cc
endif

INC = libint2.h obintlibint2.h tbintlibint2.h tform.h macros.h primpairs.h shellpairs.h gto.h

LIBSRC = $(CXXSRC) $(CSRC) $(XCSRC)
LIBOBJ = $(LIBSRC:%.c=%.$(OBJSUF))
LIBOBJ := $(LIBOBJ:%.cc=%.$(OBJSUF))

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

DISTFILES = $(TESTCSRC) $(INC) Makefile

DEPENDINCLUDE = $(INC) $(GENINC) $(SGENINC)

#############################################################

default:: $(DEPENDINCLUDE)
interface:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

TESTLIBS = $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/../wfn/LIBS.h)
libint2test: libint2test.$(OBJSUF) $(TESTLIBS) $(LIBS)
	$(LTLINK) $(CXX) $(LDFLAGS) -o libint2test $^ $(SYSLIBS) $(LTLINKBINOPTS)

libint2test.$(OBJSUF): libint2test.cc
	$(LTCOMP) $(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

targetclean::
	-rm -f libint2test

$(TESTOBJ:.$(OBJSUF)=.d) $(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) $(TESTOBJ:.$(OBJSUF)=.d)
endif

#############################################################
