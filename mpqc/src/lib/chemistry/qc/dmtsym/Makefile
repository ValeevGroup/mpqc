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

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

SGENFLAGS = -ifndef _chemistry_qc_dmtsym_sym_h

TARGET_TO_MAKE = libdmtsym
BIN_OR_LIB = LIB

CSRC = create_r.c skeleton.c syminit.c mkrpd.c mksym.c

CXXSRC = sym_cpp.cc

SGENFILE = symm.sgen
SGENSRC = symmallc.c symmasgn.c symmfree.c symmbc0.c symmrbc0.c symmsbc0.c \
          symmsnd0.c symmrcv0.c symminit.c symmiseq.c symmprnt.c symmzero.c

INC = sym_dmt.h symm_mac.h
SGENINC = $(SGENSRC:%.c=%.h)
GENINC = $(CSRC:%.c=%.gbl) $(CSRC:%.c=%.lcl)

LIBSRC = $(CSRC) $(SGENSRC) $(CXXSRC)
LIBOBJ = $(LIBSRC:%.c=%.$(OBJSUF))
LIBOBJ := $(LIBOBJ:%.cc=%.$(OBJSUF))

DISTFILES = $(CSRC) $(INC) symm.sgen Makefile LIBS.h

DEPENDINCLUDE = $(SGENINC) $(GENINC) $(INC)

#####################################################################

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

