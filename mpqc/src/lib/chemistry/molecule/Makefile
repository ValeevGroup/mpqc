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
DEFINES += -DSRCDIR=\"$(SRCDIR)\" -DSRCLIBDIR=\"$(SRCLIBDIR)\"

CXXSRC = molecule.cc molsymm.cc \
         simple.cc stre.cc bend.cc out.cc tors.cc stors.cc linip.cc linop.cc \
         energy.cc coor.cc imcoor.cc symmcoor.cc redund.cc cartcoor.cc \
         molshape.cc taylor.cc hess.cc fdhess.cc \
         molrender.cc atominfo.cc molfreq.cc formula.cc

CLASSES = Molecule \
          IntCoor SumIntCoor SetIntCoor \
          MolecularCoor IntMolecularCoor \
          MolecularEnergy  TaylorMolecularEnergy \
          VDWShape ConnollyShape \
          SimpleCo StreSimpleCo BendSimpleCo TorsSimpleCo \
          OutSimpleCo LinIPSimpleCo LinOPSimpleCo \
          AtomInfo \
          RenderedMolecule RenderedStickMolecule RenderedBallMolecule \
          RenderedMolecularSurface MolecularFormula

LIBOBJ = $(CXXSRC:%.cc=%.o)

INC = molecule.h \
      simple.h localdef.h \
      coor.h energy.h molshape.h taylor.h hess.h fdhess.h \
      molrender.h atominfo.h molfreq.h formula.h

DEPENDINCLUDE = $(INC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libmolecule

TESTSRC = moltest.cc
TESTOBJ = $(TESTSRC:%.cc=%.o)
TESTFILES = moltest.in

DISTFILES = $(CXXSRC) $(INC) Makefile $(TESTSRC) $(TESTFILES)

default:: $(DEPENDINCLUDE)

interface:: $(DEPENDINCLUDE)

LIBS = $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

LD = $(CXX)

moltest: $(TESTOBJ) $(LIBS)
	$(LD) $(LDFLAGS) -o moltest $^ $(SYSLIBS)

symmetrize: symmetrize.o $(LIBS)
	$(LD) $(LDFLAGS) -o symmetrize $^ $(SYSLIBS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

distclean::
	/bin/rm -f symmetrize moltest

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTOBJ:%.o=%.d)
endif
