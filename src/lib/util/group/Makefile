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
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

REQUIREDCXXSRC = mstate.cc \
         message.cc messimpl.cc messproc.cc \
         reduce.cc topology.cc hcube.cc \
         rnglock.cc memory.cc memmsg.cc pool.cc \
         memamsg.cc memproc.cc \
         pregtime.cc thread.cc memiter.cc memrdma.cc \
	 actmsg.cc memregion.cc

OPTIONALCXXSRC = messmpi.cc memshm.cc \
                 memmtmpi.cc

CXXSRC := $(REQUIREDCXXSRC)

ifeq ($(HAVE_SYSV_IPC),yes)
  CXXSRC += memshm.cc globcnt.cc
endif

ifeq ($(HAVE_MPI),yes)
  CXXSRC += messmpi.cc memmtmpi.cc
endif
ifeq ($(HAVE_ARMCI),yes)
  CXXSRC += memarmci.cc
endif
ifeq ($(HAVE_PTHREAD),yes)
  CXXSRC += thpthd.cc
endif

LIBOBJ = $(CXXSRC:%.cc=%.$(OBJSUF))

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libSCgroup

TESTCXXSRC = pooltest.cc messtest.cc rnglocktest.cc thrtest.cc prttest.cc actmsgtest.cc
DISTFILES = $(REQUIREDCXXSRC) $(OPTIONALCXXSRC) $(INC) \
            Makefile LIBS.h $(TESTCXXSRC)

TESTPROGS = pooltest messtest rnglocktest memtest actmsgtest

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

DEPENDINCLUDE =

default:: $(DEPENDINCLUDE)

memtest.$(OBJSUF): memtest.cc
	$(LTCOMP) $(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

prttest: prttest.$(OBJSUF) \
	libSCgroup.$(LIBSUF) \
	libSCkeyval.$(LIBSUF) \
	libSCmisc.$(LIBSUF) \
	libSCstate.$(LIBSUF) \
	libSCclass.$(LIBSUF) \
	libSCcontainer.$(LIBSUF) \
	libSCref.$(LIBSUF)
	$(LTLINK) $(LD) $(LDFLAGS) -o prttest $^ $(SYSLIBS) $(LTLINKBINOPTS)

memtest: memtest.$(OBJSUF) \
	libSCgroup.$(LIBSUF) \
	libSCkeyval.$(LIBSUF) \
	libSCmisc.$(LIBSUF) \
	libSCstate.$(LIBSUF) \
	libSCclass.$(LIBSUF) \
	libSCcontainer.$(LIBSUF) \
	libSCref.$(LIBSUF)
	$(LTLINK) $(LD) $(LDFLAGS) -o memtest $^ $(SYSLIBS) $(LTLINKBINOPTS)

pooltest: pooltest.$(OBJSUF) \
	libSCmisc.$(LIBSUF) \
	libSCgroup.$(LIBSUF) \
	libSCkeyval.$(LIBSUF) \
	libSCstate.$(LIBSUF) \
	libSCclass.$(LIBSUF) \
	libSCcontainer.$(LIBSUF) \
	libSCref.$(LIBSUF)
	$(LTLINK) $(CXX) $(CPPFLAGS) $(CXXFLAGS) -o pooltest $^ $(SYSLIBS) $(LTLINKBINOPTS)

rnglocktest: rnglocktest.$(OBJSUF) libSCmisc.$(LIBSUF) libSCgroup.$(LIBSUF)
	$(LTLINK) $(CXX) $(CPPFLAGS) $(CXXFLAGS) -o rnglocktest $^ $(SYSLIBS) $(LTLINKBINOPTS)

messtest.$(OBJSUF): messtest.cc
	$(LTCOMP) $(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

messtest: messtest.$(OBJSUF) \
	libSCgroup.$(LIBSUF) \
	libSCkeyval.$(LIBSUF) \
	libSCmisc.$(LIBSUF) \
	libSCstate.$(LIBSUF) \
	libSCclass.$(LIBSUF) \
	libSCcontainer.$(LIBSUF) \
	libSCref.$(LIBSUF)
	$(LTLINK) $(LD) $(LDFLAGS) -o messtest $^ $(SYSLIBS) $(LTLINKBINOPTS)

actmsgtest.$(OBJSUF): actmsgtest.cc
	$(LTCOMP) $(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

actmsgtest: actmsgtest.$(OBJSUF) \
	libSCgroup.$(LIBSUF) \
	libSCkeyval.$(LIBSUF) \
	libSCmisc.$(LIBSUF) \
	libSCstate.$(LIBSUF) \
	libSCclass.$(LIBSUF) \
	libSCcontainer.$(LIBSUF) \
	libSCref.$(LIBSUF)
	$(LTLINK) $(LD) $(LDFLAGS) -o actmsgtest $^ $(SYSLIBS) $(LTLINKBINOPTS)

thrtest.$(OBJSUF): thrtest.cc
	$(LTCOMP) $(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

thrtest: thrtest.$(OBJSUF) \
	libSCgroup.$(LIBSUF) \
	libSCkeyval.$(LIBSUF) \
	libSCmisc.$(LIBSUF) \
	libSCstate.$(LIBSUF) \
	libSCclass.$(LIBSUF) \
	libSCcontainer.$(LIBSUF) \
	libSCref.$(LIBSUF)
	$(LTLINK) $(LD) $(LDFLAGS) -o thrtest $^ $(SYSLIBS) $(LTLINKBINOPTS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(TESTCXXSRC:%.cc=%.d) $(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(TESTCXXSRC:%.cc=%.d) $(LIBOBJ:.$(OBJSUF)=.d)
endif
