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
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libintv2
BIN_OR_LIB = LIB

SGENFLAGS = -ifndef _chemistry_qc_intv2_libint_h

# MG is the maximum angular momentum for which we will use
# the generated build routines.
MG=none
ifeq ($(ARCH),NCUBE_V2)
  MG=2
endif
ifeq ($(ARCH),NCUBE_V3)
  MG=2
endif
ifeq ($(ARCH),SABER)
  MG=0
endif
ifeq ($(MG),none)
  MG=3
endif

TESTCSRC = inttest.cc

CSRC = bounds.c buildgc.c  comp_1e.c  comp_0e.c  normalize.c  \
       offsets.c int_fjt.c  init2e.c  atominfo.c  basis.c  comp_erep.c  \
       shiftgc.c int_print.c utils.c comp_erep23.c

CXXSRC = basisk.cc storage.cc tformv2.cc #integralv2.cc obintv2.cc tbintv2.cc

XCSRC = shellprnt.c  shellrbc0.c  shellrcv0.c  shellsbc0.c  shellsnd0.c  \
        shellasgn.c  shellzero.c

SGENFILE = atoms.sgen

SGENSRC = \
  atomsallc.c atomsbc0.c atomsfree.c \
  atomsinit.c atomsprnt.c atomsrbc0.c atomssbc0.c \
  atomsiseq.c atomsasgn.c atomssnd0.c \
  atomsrcv0.c atomszero.c

INC = bounds.h buildinter.h fjttable.h genbuild.h tformv2.h \
      storage.h int_flags.h int_libv2.h int_macros.h int_types.h inter.h \

SGENINC = $(SGENSRC:%.c=%.h)
GENINC = MG.h $(CSRC:%.c=%.gbl) $(CSRC:%.c=%.lcl)

LIBSRC = $(CSRC) $(XCSRC) $(SGENSRC) $(CXXSRC)
LIBOBJ = $(LIBSRC:%.c=%.$(OBJSUF))
LIBOBJ := $(LIBOBJ:%.cc=%.$(OBJSUF))

DISTFILES = $(TESTCSRC) $(CSRC) $(XCSRC) $(INC) atoms.sgen Makefile LIBS.h

DEPENDINCLUDE = $(INC) $(GENINC) $(SGENINC)

#############################################################

default:: $(DEPENDINCLUDE)
interface:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

inttest: inttest.o \
	$(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(CXX) $(LDFLAGS) -o inttest $^ $(SYSLIBS)

inttest.o: inttest.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

$(TESTOBJ:.o=.d) $(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTOBJ:.o=.d)
endif

#############################################################

MG.h:
	echo '#define MG $(MG)' > MG.h

#buildgc.o: buildgc.c
#	$(CC) $(CFLAGS) $(HIGHOPT_CFLAGS) $(NOALIAS_CFLAGS) -c $<
