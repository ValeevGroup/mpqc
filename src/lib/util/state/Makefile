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
LD = $(CXX)

CXXSRC = state.cc statein.cc stateout.cc state_bin.cc proxy.cc \
         state_text.cc state_file.cc translate.cc \
         tmplinst.cc

TESTPROGS = statetest

LIBOBJ = $(CXXSRC:%.cc=%.$(OBJSUF))

INC =  LIBS.h state.h statein.h stateout.h \
       translate.h \
       state_bin.h state_file.h state_text.h

DEPENDINCLUDE = $(INC) $(GENINC)

DISTFILES = $(CXXSRC) $(INC) Makefile LIBS.h

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libSCstate

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

default:: $(DEPENDINCLUDE)

interface:: $(DEPENDINCLUDE)

statetest.$(OBJSUF): statetest.cc
	$(LTCOMP) $(CXX) $(CXXFLAGS) $(DEFINES) $(INCLUDE) -DSRCDIR=\"$(SRCDIR)\" -c $<

statetest: statetest.$(OBJSUF) libSCkeyval.$(LIBSUF) \
           libSCgroup.$(LIBSUF) libSCstate.$(LIBSUF) libSCclass.$(LIBSUF) \
           libSCcontainer.$(LIBSUF) libSCref.$(LIBSUF) libSCmisc.$(LIBSUF)
	$(LTLINK) $(LD) $(LDFLAGS) -o statetest $^ $(SYSLIBS) $(LTLINKBINOPTS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

distclean::
	/bin/rm -f statetest.a.out statetest.out

statetest.d $(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) statetest.d
endif
