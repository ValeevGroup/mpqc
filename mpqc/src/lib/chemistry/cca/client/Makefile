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

default::

clean:: clean-components

CLIENT_OR_SERVER = client

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

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
include $(SRCDIR)/$(TOPDIR)/lib/cca/CCAMakefile

BIN_OR_LIB = LIB
TARGET_TO_MAKE = cca/libSCccaclient

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.$(OBJSUF)=.d): $(GEN_CCA_SRCS) $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
 include $(LIBOBJ:.$(OBJSUF)=.d)
endif
