#
# Makefile
#
# Copyright (C) 2009 Edward Valeev
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

TARGET_TO_MAKE = libSCccr12
BIN_OR_LIB = LIB

CXXSRC = ccr12.cc ccr12_info.cc ccr12_info_offset.cc \
         ccr12_info_fill_in.cc ccr12_info_fill_in_r12.cc \
         ccr12_info_routines.cc ccr12_info_drivers.cc tensor.cc mtensor.cc \
         ccsd.cc ccsd_t1.cc ccsd_t2.cc ccsd_e.cc tensorextrap.cc \
         ccsdt.cc ccsdt_t1.cc ccsdt_t2.cc ccsdt_t3.cc \
         ccsd_pt_right.cc ccsd_pt_left.cc ccsd_pt.cc \
         parenthesis2t.cc parenthesis2q.cc \
         ccsd_2t_left.cc ccsd_2t_right.cc ccsd_2q_left.cc ccsd_2q_right.cc \
         ccsdtq.cc ccsdtq_t2.cc ccsdtq_t3.cc ccsdtq_t4.cc \
         ccr12_info_drivers_r12.cc ccr12_info_intermediates.cc \
         ccsd_r12.cc ccsd_r12_e.cc ccsd_r12_t1.cc ccsd_r12_t2.cc \
         ccsdpr12.cc ccsdpr12_t1.cc ccsdpr12_t2.cc ccsdpr12_c.cc \
         ccsd_sub_r12_right.cc ccsd_sub_r12_left.cc ccsd_r12_pt_right.cc \
         lambda_ccsd_t1.cc lambda_ccsd_t2.cc lambda_ccsdpr12_t1.cc lambda_ccsdpr12_t2.cc \
         ccsd_2t_pr12_right.cc ccsd_2t_r12_left.cc \
         triples_utils.cc triples_doubles.cc triples_singles.cc triples_energy.cc

INC = 

LIBOBJ= $(CXXSRC:%.cc=%.$(OBJSUF))

DISTFILES = $(INC) Makefile LIBS.h
DEPENDINCLUDE = $(INC) $(GENINC) $(SGENINC)

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

ccr12test:: ccr12test.$(OBJSUF) $(LIBS)
	$(LTLINK) $(CXX) $(LDFLAGS) -o ccr12test $^ $(SYSLIBS) $(LTLINKBINOPTS)

ccr12test.$(OBJSUF): ccr12test.cc
	$(LTCOMP) $(CXX) $(CXXFLAGS) $(CPPFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<


ccr12test.d $(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) ccr12test.d
endif
