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
  distarray4.cc \
  distarray4_node0file.cc distarray4_memgrp.cc \
  orbitalspace.cc orbitalspace_utils.cc \
  pairiter.cc ri_basis.cc \
  transform_factory.cc transform_tbint.cc \
  transform_ijxy.cc transform_ixjy.cc transform_ikjy.cc \
  transform_iRjS.cc \
  transform_123inds.cc transform_13inds.cc transform_12inds.cc \
  compute_ixjy.cc compute_ikjy.cc compute_ijxy.cc \
  compute_iRjS.cc \
  ref.cc utils.cc \
  r12wfnworld.cc \
  twobodygrid.cc svd.cc intermeds_g12.cc \
  twoparticlecontraction.cc spin.cc \
  ebc_contribs.cc gbc_contribs.cc approxB.cc Xgen.cc FxFgen.cc \
  creator.cc approxC.cc approxCp.cc approxApp.cc \
  mp2r12_energy_util.cc compute_VXB_a.cc \
  compute_VXB_a_vbsneqobs.cc \
  gaussianfit.cc r12technology.cc twobody_intermeds.cc \
  singles_emp2.cc blas.cc fixedcoefficient.cc \
  compute_VXB_c_ansatz1.cc approxC_ansatz1.cc \
  compute_VX_GenRefansatz2.cc approxC_GenRefansatz2.cc \
  moints_runtime.cc fockbuilder.cc fockbuild_runtime.cc \
  wfnworld.cc transform_ijR.cc \
  transform_ixjy_df.cc  tbint_runtime.cc \
  r12int_eval.cc \
  coulomb.cc exchange.cc fock.cc \
  r12_amps.cc compute_amps.cc \
  mp2_pair_energies.cc \
  compute_energy_a.cc mp2r12_energy.cc bc_contrib.cc dkh.cc \
  pt2r12.cc rdm.cc
  
INC = distarray4_memgrp.h distarray4_mpiiofile.h \
distarray4_node0file.h compute_tbint_tensor.h r12int_eval.h \
container.h                     r12technology.h \
contract_tbint_tensor.h  contract_tbint_tensors_to_obtensor.h \
ref.h \
creator.h                       spin.h \
gaussianfit.h gaussianfit.timpl.h transform_123inds.h \
transform_12inds.h transform_13inds.h \
transform_factory.h \
linkage.h                       transform_ijxy.h \
mbptr12.h                       transform_ikjy.h \
orbitalspace.h                  orbitalspace_utils.h \
transform_ixjy.h \
mp2r12_energy.h                 transform_tbint.h \
mp2r12_energy_util.h            twobodygrid.h \
twobodytensorinfo.h \
pairiter.h                      twoparticlecontraction.h \
print.h                         utils.h \
print_scmat_norms.h             utils.impl.h \
r12_amps.h                      r12wfnworld.h \
distarray4.h				moints_runtime.h \
fockbuilder.h                   fockbuild_runtime.h \
transform_ijR.h \
tbint_runtime.h                 fixedcoefficient.h \
pt2r12.h			rdm.h

ifeq ($(HAVE_MPIIO),yes)
  CXXSRCS += distarray4_mpiiofile.cc
endif

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


transform_123inds.$(OBJSUF) transform_123inds.d: f77sym.h

f77sym.h: f77sym.in
	$(MKF77SYM) $< $@

distclean::
	/bin/rm -f f77sym.h

#mbptr12test.d $(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) mbptr12test.d
endif
