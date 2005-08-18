#
# Makefile
#
# Copyright (C) 2001 Edward Valeev
#
# Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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

TARGET_TO_MAKE = libSCmbptr12
BIN_OR_LIB = LIB

CXXSRCS = mbptr12.cc \
  r12ia.cc r12ia_memgrp.cc \
  r12ia_node0file.cc \
  vxb_eval_info.cc \
  compute_energy_a.cc mp2r12_energy.cc \
  multipole_ints.cc ri_basis.cc \
  moindexspace.cc pairiter.cc \
  transform_factory.cc transform_tbint.cc \
  transform_ijxy.cc transform_ixjy.cc transform_ikjy.cc \
  transform_123inds.cc transform_13inds.cc transform_12inds.cc \
  compute_ixjy.cc compute_ikjy.cc compute_ijxy.cc \
  r12int_eval.cc compute_a_gebc.cc compute_a_gebc_abs1.cc \
  coulomb.cc exchange.cc fock.cc ebc_contribs.cc gbc_contribs.cc \
  compute_a_gebc_vbs.cc compute_vxb_a_symm.cc compute_vxb_a_asymm.cc \
  dualbasis_mp2.cc r12_amps.cc compute_amps.cc \
  twobodygrid.cc svd.cc linearr12.cc intermeds_g12.cc \
  singlerefinfo.cc utils.cc compute_vxb_a_new.cc \
  twoparticlecontraction.cc spin.cc \
  mp2r12_energy_compute.cc mp2_pair_energies.cc

ifeq ($(HAVE_MPIIO),yes)
  CXXSRCS += r12ia_mpiiofile.cc
endif

LIBOBJ= $(CXXSRCS:%.cc=%.$(OBJSUF)) $(CSRCS:%.c=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d)
endif

LIBS = $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

mbptr12test:: mbptr12test.$(OBJSUF) $(LIBS)
	$(LTLINK) $(CXX) $(LDFLAGS) -o mbptr12test $^ $(SYSLIBS) $(LTLINKBINOPTS)

mbptr12test.$(OBJSUF): mbptr12test.cc
	$(LTCOMP) $(CXX) $(CXXFLAGS) $(CPPFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<


transform_123inds.$(OBJSUF) transform_123inds.d: f77sym.h

f77sym.h: f77sym.in
	$(MKF77SYM) $< $@

distclean::
	/bin/rm -f f77sym.h
