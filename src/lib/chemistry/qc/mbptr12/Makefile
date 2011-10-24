#
# Makefile
#
# Copyright (C) 2001 Edward Valeev
#
# Author: Edward Valeev <evaleev@vt.edu>
# Maintainer: EV
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

TARGET_TO_MAKE = libSCmbptr12
BIN_OR_LIB = LIB

CXXSRCS = mbptr12.cc \
  ri_basis.cc \
  ref.cc \
  r12wfnworld.cc \
  twobodygrid.cc intermeds_g12.cc \
  twoparticlecontraction.cc \
  ebc_contribs.cc gbc_contribs.cc approxB.cc Xgen.cc FxFgen.cc \
  creator.cc approxC.cc approxCp.cc approxApp.cc \
  mp2r12_energy_util.cc compute_VXB_a.cc \
  compute_VXB_a_vbsneqobs.cc \
  gaussianfit.cc r12technology.cc twobody_intermeds.cc \
  singles_emp2.cc fixedcoefficient.cc \
  compute_VXB_c_ansatz1.cc approxC_ansatz1.cc \
  compute_VX_GenRefansatz2.cc approxC_GenRefansatz2.cc \
  r12int_eval.cc \
  coulomb.cc exchange.cc fock.cc \
  r12_amps.cc compute_amps.cc \
  mp2_pair_energies.cc \
  compute_energy_a.cc mp2r12_energy.cc bc_contrib.cc dkh.cc \
  pt2r12.cc pt2r12_spinorbital.cc
  
INC = compute_tbint_tensor.h r12int_eval.h \
container.h                     r12technology.h \
contract_tbint_tensor.h  contract_tbint_tensors_to_obtensor.h \
ref.h \
creator.h \
gaussianfit.h gaussianfit.timpl.h \
linkage.h mbptr12.h \
mp2r12_energy.h mp2r12_energy_util.h            twobodygrid.h \
twobodytensorinfo.h \
pairiter.h                      twoparticlecontraction.h \
print.h \
r12_amps.h                      r12wfnworld.h \
fixedcoefficient.h              pt2r12.h

ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_CCA),yes)
  DEFINES += -DHAVE_CHEMISTRY_CCA
endif

LIBOBJ= $(CXXSRCS:%.cc=%.$(OBJSUF)) $(CSRCS:%.c=%.$(OBJSUF))

DISTFILES = $(INC) Makefile LIBS.h
DEPENDINCLUDE = $(INC) $(GENINC) $(SGENINC)

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

default::

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

mbptr12test:: mbptr12test.$(OBJSUF) $(LIBS)
	$(LTLINK) $(CXX) $(LDFLAGS) -o mbptr12test $^ $(SYSLIBS) $(LTLINKBINOPTS)

mbptr12test.$(OBJSUF): mbptr12test.cc
	$(LTCOMP) $(CXX) $(CXXFLAGS) $(CPPFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

#mbptr12test.d $(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) mbptr12test.d
endif
