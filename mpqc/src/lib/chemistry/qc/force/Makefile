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
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libforce
BIN_OR_LIB = LIB

DEFINES += -DTIMING

HISTORICALWEBSRC = spcsscf.web mpcsscf.web mposscf.web
WEBSRC = mpcsscf.web mposscf.web
.PRECIOUS: $(WEBSRC:%.web=%.c)
CXXSRC = read.cc

GENSRC = $(WEBSRC:%.web=%.c)
SRCS = $(GENSRC)

INC = libforce.h read.h
GENINC = $(SRCS:%.c=%.gbl)

DEPENDINCLUDE = $(INC) $(GENINC)

LIBOBJ= $(SRCS:%.c=%.$(OBJSUF)) $(CXXSRC:%.cc=%.$(OBJSUF))

################################################################3

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

ifneq ($(SRCDIR),.)
%.c: %.web
	cp $< .
	$(CTANGLE) $(notdir $<)
	rm -f $(notdir $<)

%.tex: %.web
	cp $< .
	$(CWEAVE) $(notdir $<)
	rm -f $(notdir $<)
else
%.c: %.web
	$(CTANGLE) $<

%.tex: %.web
	$(CWEAVE) $<
endif
