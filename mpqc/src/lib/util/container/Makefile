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

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

CXXSRC = array.cc avl.cc bitarray.cc carray.cc

TESTCXX = avltest.cc artest.cc

TESTPROGS = avltest artest

LIBOBJ = $(CXXSRC:%.cc=%.$(OBJSUF))

INC = array.h artem.h ar2tem.h ssarmac.h ssar2mac.h \
      eavlmmap.h avlmap.h avlset.h compare.h bitarray.h

GENINC = armac.h ar2mac.h ssarmac.h ssar2mac.h

DEPENDINCLUDE = $(INC) $(GENINC)

DISTFILES = Makefile $(INC) LIBS.h $(REALCXXSRC) $(TESTCXX)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libSCcontainer

default:: $(DEPENDINCLUDE)

interface:: $(DEPENDINCLUDE)

armac.h: artem.h $(TMPL2MACDEP)
	$(TMPL2MAC) $< > armac.h

ar2mac.h: ar2tem.h $(TMPL2MACDEP)
	$(TMPL2MAC) $< > ar2mac.h

ssarmac.h: ssartem.h $(TMPL2MACDEP)
	$(TMPL2MAC) $< > ssarmac.h

ssar2mac.h: ssar2tem.h $(TMPL2MACDEP)
	$(TMPL2MAC) $< > ssar2mac.h

avltest: avltest.$(OBJSUF) libSCmisc.$(LIBSUF) libSCkeyval.$(LIBSUF) libSCstate.$(LIBSUF) libSCclass.$(LIBSUF) libSCref.$(LIBSUF) libSCcontainer.$(LIBSUF)
	$(LTLINK) $(CXX) $(LDFLAGS) -o avltest $^ $(SYSLIBS) $(LTLINKBINOPTS)

artest: artest.$(OBJSUF) libSCcontainer.$(LIBSUF) libSCmisc.$(LIBSUF) libSCkeyval.$(LIBSUF) libSCstate.$(LIBSUF) libSCclass.$(LIBSUF) libSCref.$(LIBSUF) libSCcontainer.$(LIBSUF)
	$(LTLINK) $(CXX) $(LDFLAGS) -o artest $^ $(SYSLIBS) $(LTLINKBINOPTS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

miscclean::
	-rm -f foo.h

$(TESTCXX:%.cc=%.d) $(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(TESTCXX:%.cc=%.d) $(LIBOBJ:.$(OBJSUF)=.d)
endif
