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
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

REQUIREDCXXSRC = mstate.cc \
         message.cc messagei.cc messproc.cc messint.cc \
         reduce.cc topology.cc hcube.cc picl.cc \
         rnglock.cc memory.cc memmsg.cc pool.cc memoryi.cc \
         globcnt.cc memamsg.cc memproc.cc memmid.cc \
         pregtime.cc

OPTIONALCXXSRC = messshm.cc messpgon.cc messpvm.cc messmpi.cc memshm.cc \
                 mempgon.cc memipgon.cc memmpl.cc

CXXSRC := $(REQUIREDCXXSRC)

ifeq ($(HAVE_SYSV_IPC),yes)
  CXXSRC += messshm.cc memshm.cc
endif
ifeq ($(HAVE_NX),yes)
  CXXSRC += messpgon.cc memipgon.cc
ifeq ($(HAVE_HRECV),yes)
  CXXSRC += mempgon.cc
endif
endif
ifeq ($(HAVE_PVM),yes)
  CXXSRC += messpvm.cc
endif
ifeq ($(HAVE_MPI),yes)
  CXXSRC += messmpi.cc memmpi.cc
endif
ifeq ($(HAVE_MPL),yes)
  CXXSRC += memmpl.cc
endif

INC = pool.h message.h mstate.h topology.h picl.h \
      messshm.h messpvm.h messmpi.h \
      rnglock.h memory.h memmsg.h memshm.h pool.h memamsg.h memproc.h \
      memmpl.h memmpi.h pregtime.h

LIBOBJ = $(CXXSRC:%.cc=%.o)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libgroup

TESTCXXSRC = pooltest.cc messtest.cc rnglocktest.cc
DISTFILES = $(REQUIREDCXXSRC) $(OPTIONALCXXSRC) $(INC) \
            Makefile LIBS.h $(TESTCXXSRC)

DEPENDINCLUDE =

default:: $(DEPENDINCLUDE)

memtest.o: memtest.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

memtest: memtest.o \
	libgroup.a \
	libkeyval.a \
	libmisc.a \
	libstate.a \
	libclass.a \
	libcontainer.a \
	libref.a
	$(LD) $(LDFLAGS) -o memtest $^ $(SYSLIBS)

pooltest: pooltest.o libmisc.a libgroup.a libkeyval.a libstate.a \
          libclass.a libcontainer.a libref.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o pooltest $^ $(SYSLIBS)

rnglocktest: rnglocktest.o libgroup.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o rnglocktest $^ $(SYSLIBS)

messtest.o: messtest.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

messtest: messtest.o \
	libgroup.a \
	libkeyval.a \
	libmisc.a \
	libstate.a \
	libclass.a \
	libcontainer.a \
	libref.a
	$(LD) $(LDFLAGS) -o messtest $^ $(SYSLIBS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(TESTCXXSRC:%.cc=%.d) $(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(TESTCXXSRC:%.cc=%.d) $(LIBOBJ:.o=.d)
endif

