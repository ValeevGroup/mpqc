#
# Makefile
#
# Copyright (C) 1996 Limit Point Systems, Inc.
#
# Author: Edward Seidl <seidl@janed.com>
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

DEFINES += -DSRCDIR=\"$(SRCDIR)\"

TARGET_TO_MAKE = libSCbasis
BIN_OR_LIB = LIB

CXXSRC = gaussbas.cc gaussshell.cc files.cc \
         gaussbaval.cc gaussshval.cc cartiter.cc \
         transform.cc shellrot.cc integral.cc obint.cc petite.cc aotoso.cc \
         symmint.cc tbint.cc dercent.cc \
         sobasis.cc sointegral.cc extent.cc orthog.cc \
         gpetite.cc distshpair.cc intparams.cc \
         uncontract.cc split.cc lselect.cc operator.cc fjt.cc \
         distsh.cc operator.cc union.cc

INC = basis.h gaussbas.h guassshell.h files.h \
      cartiter.h transform.h shellrot.h \
      integral.h obint.h petite.h symmint.h tbint.h dercent.h \
      sobasis.h sointegral.h \
      orthog.h gpetite.h distshpair.h intparams.h intdescr.h \
      uncontract.h split.h lselect.h operator.h fjt.h \
      distsh.h inttraits.h operator.h union.h \
      obintfactory.h

TESTSRC = btest.cc tbint_batch_test.cc
TESTOBJ = $(TESTSRC:%.cc=%.$(OBJSUF))
TESTFILES = btest.kv tbint_batch_test.in

TESTPROGS = btest tbint_batch_test

DISTFILES = $(CXXSRC) $(INC) Makefile $(TESTSRC) $(TESTFILES)

LIBOBJ= $(CXXSRC:%.cc=%.$(OBJSUF))

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

default:: $(DEPENDINCLUDE)

btest: $(TESTOBJ) libSCintv3.$(LIBSUF) libSCoint3.$(LIBSUF) \
                  $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LTLINK) $(LD) $(LDFLAGS) -o btest $^ $(SYSLIBS) $(LTLINKBINOPTS)

tbint_batch_test:: tbint_batch_test.$(OBJSUF) libSCscf.$(LIBSUF) libSCwfn.$(LIBSUF) libSCintv3.$(LIBSUF) libSCoint3.$(LIBSUF) \
                  $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LTLINK) $(CXX) $(LDFLAGS) -o tbint_batch_test $^ $(SYSLIBS) $(LTLINKBINOPTS)

tbint_batch_test.$(OBJSUF): tbint_batch_test.cc
	$(LTCOMP) $(CXX) $(CXXFLAGS) $(CPPFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

clean::
	rm -f btest btest.$(OBJSUF)ut

$(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) $(TESTOBJ:%.$(OBJSUF)=%.d)
endif

