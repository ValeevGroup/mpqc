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
include $(TOPDIR)/lib/Makedirlist

SUBDIRS = basis oint3 intv3 wfn scf dft mbpt
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_CC),yes)
SUBDIRS := $(SUBDIRS) cc
endif
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_CINTS),yes)
SUBDIRS := $(SUBDIRS) cints
endif
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_LIBINT2),yes)
SUBDIRS := $(SUBDIRS) libint2
endif
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_PSI),yes)
SUBDIRS := $(SUBDIRS) psi
endif
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_MBPTR12),yes)
SUBDIRS := $(SUBDIRS) mbptr12
endif
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_INTCCA),yes)
SUBDIRS := $(SUBDIRS) intcca
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalSubDirs

