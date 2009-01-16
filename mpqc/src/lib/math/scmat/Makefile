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
DEFINES += -DSRCDIR=\"$(SRCDIR)\"

CLASS2TEXFLAGS = -clssection subsection -clssubsection subsubsection

default:: $(DEPENDINCLUDE)

CXXSRC = dim.cc \
         abstract.cc   matrix.cc     matrix_i.cc \
         block.cc      blkiter.cc    elemop.cc     result.cc \
         local.cc      localrect.cc  localdiag.cc  localsymm.cc  localvect.cc \
         disthql.cc \
         dist.cc       distrect.cc   distdiag.cc   distsymm.cc   distvect.cc \
         repl.cc       replrect.cc   repldiag.cc   replsymm.cc   replvect.cc \
         blocked.cc \
         blockedvect.cc blockedrect.cc blockeddiag.cc blockedsymm.cc \
         matrix3.cc    vector3.cc    vector3_i.cc \
         util.cc       stateutil.cc

TESTSRC = disttest.cc localtest.cc repltest.cc blockedtest.cc

TESTPROGS = disttest localtest repltest blockedtest

CSRC = cmatrix.c svd.c

FSRC = pdsteqr.f

LIBOBJ = $(CXXSRC:%.cc=%.$(OBJSUF)) $(CSRC:%.c=%.$(OBJSUF)) $(FSRC:%.f=%.$(OBJSUF)) 

INC =  cmatrix.h matrix_i.h cmatrix.h matrix.h abstract.h blkiter.h local.h \
       result.h vector3.h vector3_i.h block.h matrix3.h elemop.h repl.h \
       blocked.h disthql.h dim.h

GENINC =

DEPENDINCLUDE = $(INC) $(GENINC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libSCscmat

TESTOBJ = matrixtest.$(OBJSUF)

DISTFILES = $(CXXSRC) $(INC) scmat.h $(WEBSRC) Makefile LIBS.h \
            $(TESTSRC) $(TESTFILES)

LIBS = $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

matrixtest:
	@echo The target must be a specialization such as localtest.

localtest: localtest.$(OBJSUF) $(TESTOBJ) $(LIBS)
	$(LTLINK) $(LD) $(LDFLAGS) -o localtest $^ $(SYSLIBS) $(LTLINKBINOPTS)

repltest: repltest.$(OBJSUF) $(TESTOBJ) $(LIBS)
	$(LTLINK) $(LD) $(LDFLAGS) -o repltest $^ $(SYSLIBS) $(LTLINKBINOPTS)

disttest: disttest.$(OBJSUF) $(TESTOBJ) $(LIBS)
	$(LTLINK) $(LD) $(LDFLAGS) -o disttest $^ $(SYSLIBS) $(LTLINKBINOPTS)

blockedtest: blockedtest.$(OBJSUF) $(TESTOBJ) $(LIBS)
	$(LTLINK) $(LD) $(LDFLAGS) -o blockedtest $^ $(SYSLIBS) $(LTLINKBINOPTS)

repltest.$(OBJSUF): repltest.cc
	$(LTCOMP) $(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $< -o $@

disttest.$(OBJSUF): disttest.cc
	$(LTCOMP) $(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $< -o $@

TAGS: $(CXXSRC) $(CSRC)
	etags $^

disthql.$(OBJSUF) disthql.d: f77sym.h

f77sym.h: f77sym.in
	$(MKF77SYM) $< $@

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

distclean::
	/bin/rm -f f77sym.h

$(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(CXXSRC:.cc=.d) $(CSRC:.c=.d) $(TESTSRC:%.cc=%.d)
endif

