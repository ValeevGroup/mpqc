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

TARGET_TO_MAKE = libmbptr12
BIN_OR_LIB = LIB

CXXSRCS = mbptr12.cc distsh.cc trans12_grt.cc \
	  trans123_r12a_abs.cc r12ia.cc r12ia_memgrp.cc \
	  r12ia_node0file.cc vxb_eval.cc vxb_eval_info.cc \
	  vxb_eval_sbs_a.cc compute_sbs_a.cc \
	  vxb_eval_abs_a.cc compute_abs_a.cc \
	  compute_energy_a.cc mp2r12_energy.cc \
	  multipole_ints.cc

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

mbptr12test:: mbptr12test.$(OBJSUF) $(LIBS) libdft.$(LIBSUF)
	$(LTLINK) $(CXX) $(LDFLAGS) -o mbptr12test $^ $(SYSLIBS) $(LTLINKBINOPTS)

mbptr12test.$(OBJSUF): mbptr12test.cc
	$(LTCOMP) $(CXX) $(CXXFLAGS) $(CPPFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<
