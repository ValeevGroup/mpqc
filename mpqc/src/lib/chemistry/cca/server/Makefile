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
include $(TOPDIR)/lib/Makedirlist

default::

clean:: clean-components

CLIENT_OR_SERVER = server

BABEL_PACKAGES = MPQC

BABEL_CLASSES = \
  MPQC.CCAException \
  MPQC.Chemistry_MoleculeViewer \
  MPQC.Chemistry_QC_Model \
  MPQC.Chemistry_QC_ModelFactory \
  MPQC.Physics_Units \
  MPQC.GaussianBasis_Shell \
  MPQC.GaussianBasis_Atomic \
  MPQC.GaussianBasis_Molecular \
  MPQC.IntV3EvaluatorFactory \
  MPQC.IntegralEvaluator1 \
  MPQC.IntegralEvaluator2 \
  MPQC.IntegralEvaluator3 \
  MPQC.IntegralEvaluator4 \
  MPQC.ChemistryOpt_CoordinateModel \
  MPQC.ComponentClassDescription \
  MPQC.ComponentFactory \
  MPQC.SimpleDriver

ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_CINTS),yes)
  BABEL_CLASSES += MPQC.CintsEvaluatorFactory
endif

EXTRA_CXX_SRCS = socket.cc except.cc basis_cca_to_sc.cc reorder_engine.cc \
                 buffer_size.cc
EXTRA_INCLUDES = socket.h except.h ccaiter.h \
                 integral_evaluator.h reorder_engine.h buffer_size.h

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
include $(SRCDIR)/$(TOPDIR)/lib/cca/CCAMakefile

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libSCccaserver

DEFINES += -DSRCDIR=\"$(SRCDIR)\"
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_CINTS),yes)
  DEFINES += -DHAVE_CINTS
endif

LTLINKLIBOPTS += -L$(libdir) -R$(libdir)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.$(OBJSUF)=.d): $(GEN_CCA_SRCS) $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
 include $(LIBOBJ:.$(OBJSUF)=.d)
endif
