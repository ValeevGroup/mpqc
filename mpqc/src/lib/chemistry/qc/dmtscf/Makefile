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
  SRCDIR=.
endif

SGENFLAGS = -ifndef _chemistry_qc_dmtscf_scf_h

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libdmtscf
BIN_OR_LIB = LIB

CSRC = scf_bnd.c scf_core.c scf_diis.c scf_en.c \
       scf_gmat.c scf_iter.c scf_loopg.c \
       scf_mkden.c scf_mkgd.c scf_oeis.c scf_orth.c \
       scf_proj.c scf_vect.c scf_ex.c scf_mkgdlb.c scf_mkk.c scf_mkj.c \
       scf_loopj.c scf_loopk.c

CXXSRC = scf_kinp.cc

SGENFILE = scf.sgen

SGENSRC = scfasgn.c scfbc0.c scfinit.c \
          scfiseq.c scfprnt.c scfrbc0.c scfsbc0.c scfsnd0.c scfrcv0.c \
          scfallc.c scffree.c scfzero.c
GENSRC = $(SGENSRC)

INC = scf_dmt.h

SGENINC = $(GENSRC:%.c=%.h)
GENINC = $(CSRC:%.c=%.gbl) $(CSRC:%.c=%.lcl)

LIBOBJ = $(CSRC:%.c=%.o) $(CXXSRC:%.cc=%.o) $(GENSRC:%.c=%.o)

DEPENDINCLUDE = $(INC) $(GENINC) $(SGENINC)

DISTFILES = $(CSRC) $(INC) scf.sgen Makefile LIBS.h

###################################################################

default:: $(DEPENDINCLUDE)
interface:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

dmtscftest: dmtscftest.o libdmtscf.a \
        $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(CXX) $(LDFLAGS) -o $@ $^ $(SYSLIBS)

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) dmtscftest.d
endif

###################################################################

# The SGI cc compiler cannot handle scf_mkgd.c.
# The resulting code will sometimes produce a wrong answer.
#ifeq ($(ARCH),SGI)
#scf_mkgd.o: scf_mkgd.c
#	gcc -O -g $(INCLUDE) -c $< -o scf_mkgd.o
#scf_ex.o: scf_ex.c
#	gcc -g $(INCLUDE) -c $< -o scf_ex.o
#endif
