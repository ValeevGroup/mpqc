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

CXXSRC = avl.cc bitarray.cc carray.cc

TESTCXX = avltest.cc

TESTPROGS = avltest

LIBOBJ = $(CXXSRC:%.cc=%.$(OBJSUF))

INC = eavlmmap.h avlmap.h avlset.h compare.h bitarray.h

DEPENDINCLUDE = $(INC)

DISTFILES = Makefile $(INC) LIBS.h $(REALCXXSRC) $(TESTCXX)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libSCcontainer

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

default:: $(DEPENDINCLUDE)

interface:: $(DEPENDINCLUDE)

avltest: avltest.$(OBJSUF) libSCmisc.$(LIBSUF) libSCkeyval.$(LIBSUF) libSCstate.$(LIBSUF) libSCclass.$(LIBSUF) libSCref.$(LIBSUF) libSCcontainer.$(LIBSUF)
	$(LTLINK) $(CXX) $(LDFLAGS) -o avltest $^ $(SYSLIBS) $(LTLINKBINOPTS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

miscclean::
	-rm -f foo.h

$(TESTCXX:%.cc=%.d) $(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(TESTCXX:%.cc=%.d) $(LIBOBJ:.$(OBJSUF)=.d)
endif
