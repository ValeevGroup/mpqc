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

default::

clean:: clean-components

BABEL_PACKAGES = MPQC

BABEL_CLASSES = \
  MPQC.Chemistry_MoleculeViewer \
  MPQC.Chemistry_QC_Model \
  MPQC.Chemistry_QC_ModelFactory \
  MPQC.Physics_Units \
  MPQC.GaussianBasis_Shell \
  MPQC.GaussianBasis_Atomic \
  MPQC.GaussianBasis_Molecular \
  MPQC.IntV3EvaluatorFactory \
  MPQC.CintsEvaluatorFactory \
  MPQC.IntegralEvaluator1 \
  MPQC.IntegralEvaluator2 \
  MPQC.IntegralEvaluator3 \
  MPQC.IntegralEvaluator4 \
  MPQC.ChemistryOpt_CoordinateModel \
  MPQC.ComponentClassDescription \
  MPQC.ComponentFactory \
  MPQC.SimpleDriver

EXTRA_CXX_SRCS = socket.cc except.cc basis_cca_to_sc.cc integral_evaluator.cc
EXTRA_INCLUDES = socket.h except.h ccaiter.h integral_evaluator.h

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
DEFINES += -DSRCDIR=\"$(SRCDIR)\"

SIDL_FILES=$(SRCTOPDIR)/lib/cca/sidl/chemistry-mpqc.sidl
REPO = $(TOPDIR)/lib/cca/repo
BABEL_ARGS = -E --suppress-timestamp -R$(CCA_SPEC_BABEL_SHARE)/xml -R$(CCA_CHEM_REPO) -R$(REPO) --server=C++

## so we don't have to include babel.make
PACKAGE_BASES := $(shell echo $(BABEL_PACKAGES) | sed s/\\./\\_/g )
STUBHDRS = $(PACKAGE_BASES:%=%.hh)
IORHDRS = $(PACKAGE_BASES:%=%_IOR.h)

ENUM_BASES := $(shell echo $(BABEL_ENUMS) | sed s/\\./\\_/g )
STUBHDRS += $(ENUM_BASES:%=%.hh)
IORHDRS  += $(ENUM_BASES:%=%_IOR.h)
IORSRCS  = $(ENUM_BASES:%=%_IOR.c)

CLASS_BASES := $(shell echo $(BABEL_CLASSES) | sed s/\\./\\_/g )
IORHDRS += $(CLASS_BASES:%=%_IOR.h)
IORSRCS += $(CLASS_BASES:%=%_IOR.c)
STUBHDRS += $(CLASS_BASES:%=%.hh)
STUBSRCS = $(CLASS_BASES:%=%.cc)
IMPLHDRS = $(CLASS_BASES:%=%_Impl.hh)
IMPLSRCS = $(CLASS_BASES:%=%_Impl.cc)
SKELSRCS = $(CLASS_BASES:%=%_Skel.cc)

GEN_SRCS = $(IORHDRS) $(IORSRCS) $(STUBHDRS) $(STUBSRCS) \
           $(IMPLHDRS) $(IMPLSRCS) $(SKELSRCS)
CCA_CXX_SRCS = $(IMPLSRCS) $(SKELSRCS) $(STUBSRCS) $(EXTRA_CXX_SRCS)
CCA_CC_SRCS = $(IORSRCS)
CCA_OBJS = $(CCA_CXX_SRCS:%.cc=%.$(OBJSUF)) $(CCA_CC_SRCS:%.c=%.$(OBJSUF))
CCA_INCLUDES = $(STUBHDRS) $(IORHDRS) $(EXTRA_INCLUDES)
###

BIN_OR_LIB = LIB
TARGET_TO_MAKE = cca/libSCchemistry

LIBOBJ = $(CCA_OBJS)
CPPFLAGS += -I./ -I$(SRCDIR) -I$(CCA_SPEC_CLASSIC_INCLUDE)
ifeq ($(HAVE_LIBINT),yes)
  DEFINES += -DHAVE_CINTS
endif
ifeq ($(INTV3_ORDER),yes)
  DEFINES += -DINTV3_ORDER
endif
tmp_libdir := $(libdir)
SCLIBS = $(shell echo `ls ../../../../lib/*.$(LIBSUF)`)
SCLINK = $(SCLIBS:lib%.$(LIBSUF)=-l%)
LTLINKLIBOPTS += -L$(CCA_CHEM_LIB) -R$(CCA_CHEM_LIB) -lccachem_cxx_client -lccachem_cxx_server -L$(tmp_libdir) -R$(tmp_libdir) $(SCLINK) $(SYSLIBS)
libdir = $(tmp_libdir)/cca

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d)
endif
