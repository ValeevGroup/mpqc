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

TARGET_TO_MAKE=mpqc
# This is required to give libtool the right rpath arguments.
# Without the right rpath arguments, libtool will refuse to install
# the module.
libdir:=$(libdir)/cca

default:: $(TOPDIR)/lib/cca/$(TARGET_TO_MAKE).la

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

SRC=dummy.cc
LIBOBJ=$(SRC:%.cc=%.$(OBJSUF))

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

$(TOPDIR)/lib/cca/$(TARGET_TO_MAKE).la: $(LIBOBJ) $(LIBS)
	$(LTLINK) $(CXX) -module -o $@ $^ $(LTLINKLIBOPTS) $(SYSLIBS)

install:: $(TOPDIR)/lib/cca/$(TARGET_TO_MAKE).$(LIBSUF)
	$(INSTALL) $(INSTALLDIROPT) $(installroot)$(libdir)
	$(LTINST) $(INSTALL) -module $(INSTALLLIBOPT) $< $(installroot)$(libdir)
