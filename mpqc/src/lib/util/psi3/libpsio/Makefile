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

CXXSRC = init.cc done.cc open.cc close.cc open_check.cc error.cc \
rw.cc read.cc write.cc \
tocwrite.cc tocread.cc tocscan.cc toclast.cc tocdel.cc tocclean.cc tocprint.cc \
get_filename.cc get_numvols.cc get_volpath.cc volseek.cc get_address.cc get_global_address.cc \
write_entry.cc read_entry.cc write_block.cc read_block.cc

LIBOBJ = $(CXXSRC:%.cc=%.$(OBJSUF)) $(GENCXXSRC:%.cc=%.$(OBJSUF))

INC = psio.h psio.gbl

DEPENDINCLUDE = $(INC) $(GENINC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libSCpsio

TESTCXXSRC = psiotest.cc
DISTFILES = $(CXXSRC) $(INC) Makefile LIBS.h $(TESTCXXSRC)

TESTPROGS = psiotest

default:: $(DEPENDINCLUDE)

interface:: $(DEPENDINCLUDE)

# a test program
LD = $(CXX)
psiotest: psiotest.$(OBJSUF) libSCpsio.$(LIBSUF) libSCclass.$(LIBSUF) libSCcontainer.$(LIBSUF) libSCref.$(LIBSUF) libSCmisc.$(LIBSUF) libSCstate.$(LIBSUF) libSCkeyval.$(LIBSUF)
	$(LTLINK) $(LD) $(LDFLAGS) -o psiotest $^ $(SYSLIBS) $(LTLINKBINOPTS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(TESTCXXSRC:%.cc=%.d) $(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(TESTCXXSRC:%.cc=%.d) $(LIBOBJ:.$(OBJSUF)=.d)
endif
