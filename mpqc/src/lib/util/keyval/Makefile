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
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
LD = $(CXX)

CXXSRC = ipv2.cc ipv2_alloc.cc ipv2_cwk.cc ipv2_data.cc ipv2_error.cc \
	 ipv2_karray.cc ipv2_print.cc ipv2_read.cc \
	 keyval.cc keyvalipv2.cc keyvalval.cc keyvalkey.cc keyvalass.cc \
         keyvalstr.cc keyvalagg.cc keyvalpre.cc

GENCXXSRC =ipv2_parse.cc ipv2_scan.cc

LIBOBJ = $(CXXSRC:%.cc=%.o) $(GENCXXSRC:%.cc=%.o)
GENSRC = $(GENCXXSRC)

INC = ipv2.h keyval.h keyvaldefs.h ipv2_scan.h
GENINC = ipv2_parse.h


DEPENDINCLUDE = $(INC) $(GENINC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libkeyval

DISTFILES = Makefile $(CXXSRC) $(INC) ipv2_parse.yy ipv2_scan.ll LIBS.h

default:: $(DEPENDINCLUDE)

keyvaltest: keyvaltest.o libkeyval.$(LIBSUF) libcontainer.$(LIBSUF) \
            libmisc.$(LIBSUF) libgroup.$(LIBSUF) \
            libclass.$(LIBSUF) libstate.$(LIBSUF) \
            libref.$(LIBSUF)
	$(LD) $(LDFLAGS) -o keyvaltest $^ $(SYSLIBS)

keyvaltest.o: keyvaltest.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

distclean::
	/bin/rm -f keyvaltest

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
OTHEROBJ = keyvaltest.o
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(OTHEROBJ:.o=.d)
endif

#### yacc and lex ####
# (only works with bison and flex)
ipv2_parse.cc ipv2_parse.h: ipv2_parse.yy
	$(BISON) -v -d -o ipv2_parse.tmp $^
	sed "s/^int yyparse.*;$$//" < ipv2_parse.tmp > ipv2_parse.cc
	/bin/mv ipv2_parse.tmp.h ipv2_parse.h
	-@rm -f ipv2_parse.tmp
	-@rm -f ipv2_parse.output

ipv2_scan.cc: ipv2_scan.ll
	$(FLEX) -L -t $^ | grep -v "extern FILE .yyin" \
	                | grep -v "static int yy_get_next_buffer.*;" \
	                | sed "s/static int yy_get_next_buffer/int yy_get_next_buffer/" \
	                | grep -v "static void yyunput.*;" \
	                | sed "s/static void yyunput/void yyunput/" \
	                | grep -v "static int yyinput.*;" \
	                | sed "s/static int yyinput/int yyinput/" \
	                > ipv2_scan.cc

