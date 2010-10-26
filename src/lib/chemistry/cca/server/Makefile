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
include $(TOPDIR)/lib/Makedirlist

ifeq ($(INTV3_ORDER),yes)
  DEFINES += -DINTV3_ORDER
endif

DEFINES += -DSRCDIR=\"$(SRCDIR)\"
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_CINTS),yes)
  DEFINES += -DHAVE_CINTS
endif
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_LIBINT2),yes)
  DEFINES += -DHAVE_LIBINT2
endif

include babel.make

TARGET_TO_MAKE = libSCccaserver
BIN_OR_LIB = LIB

CSRC = $(IORSRCS)

CXXSRC = $(IMPLSRCS) $(SKELSRCS) $(STUBSRCS) \
  basis_cca_to_sc.cc buffer_size.cc except.cc reorder_engine.cc socket.cc

INC = $(STUBHDRS) 

CPPFLAGS += -I./

LTLINKLIBOPTS += \
  -L$(libdir) -R$(libdir) \
  -L$(BABEL_LIB) -R$(BABEL_LIB) -lsidl \
  -L$(CCAFE_LIB) -R$(CCAFE_LIB) $(CCAFE_CORE) \
  -L$(CCA_CHEM_LIB) -R$(CCA_CHEM_LIB) -lccachem

LIBSRC = $(CXXSRC) $(CSRC)
LIBOBJ = $(LIBSRC:%.c=%.$(OBJSUF))
LIBOBJ := $(LIBOBJ:%.cc=%.$(OBJSUF))
LIBOBJ := $(LIBOBJ:%.cxx=%.$(OBJSUF))

DISTFILES = $(INC) Makefile

DEPENDINCLUDE = $(INC)

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

#############################################################

default:: $(DEPENDINCLUDE)
interface:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

babel.make: $(SRCDIR)/$(TOPDIR)/lib/cca/sidl/chemistry-mpqc.sidl
	cd ..; $(MAKE) babel-mpqc-cxx

targetclean::

distclean::
	-rm -f $(IORHDRS) $(IORSRCS) $(SKELSRCS) $(STUBHDRS) $(STUBSRCS)
	-rm -f babel*

$(TESTOBJ:.$(OBJSUF)=.d) $(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) $(TESTOBJ:.$(OBJSUF)=.d)
endif

