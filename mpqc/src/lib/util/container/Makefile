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

SETIMPL = AVL
MAPIMPL = RAVL

GENCXXSRC = pixvpMap.cc pixvp$(MAPIMPL)Map.cc \
            voidptrSet.cc voidptr$(SETIMPL)Set.cc \
            pixintMap.cc pixint$(MAPIMPL)Map.cc \
            intpixMap.cc intpix$(MAPIMPL)Map.cc

REALCXXSRC = array.cc set.cc voidptr.cc avl.cc error.cc

GENSRC = $(GENCXXSRC)

CXXSRC = $(GENCXXSRC) $(REALCXXSRC)
TESTCXX = settest.cc settestt.cc

LIBOBJ = $(CXXSRC:%.cc=%.o)

INC = pixmap.h pixvpdefs.h ptrset.h voidptr.h voidptrdefs.h \
      set.h array.h artem.h ar2tem.h ssarmac.h ssar2mac.h \
      pixintdefs.h avl.h refset.h

GENINC = $(GENCXXSRC:%.cc=%.h) setmacr.h asetmacr.h \
         armac.h ar2mac.h ssarmac.h ssar2mac.h

DEPENDINCLUDE = $(INC) $(GENINC)

DISTFILES = Makefile $(INC) LIBS.h $(REALCXXSRC) $(TESTCXX)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libcontainer

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

setmacr.h: settmpl.h $(TMPL2MACDEP)
	$(TMPL2MAC) $< > setmacr.h

asetmacr.h: asettmpl.h $(TMPL2MACDEP)
	$(TMPL2MAC) $< > asetmacr.h

settest: settest.o libcontainer.$(LIBSUF) libmisc.$(LIBSUF) libref.$(LIBSUF)
	$(CXX) $(LDFLAGS) -o settest $^

artest: artest.o libcontainer.$(LIBSUF)
	$(CXX) $(LDFLAGS) -o artest $^

settestt: settestt.o libcontainer.$(LIBSUF)
	$(CXX) $(LDFLAGS) -o settestt $^

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

miscclean::
	-rm -f foo.h

$(TESTCXX:%.cc=%.d) $(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

pixintMap.cc pixintMap.h:
	$(GENCLASS) -2 Pix val int val Map pixint
	$(PREPEND_HEADER) pixintMap.h "#include <Pix.h>"
	$(PREPEND_HEADER) pixintMap.h "#include <util/container/pixintdefs.h>"
	/bin/mv pixintMap.h tmp.pixintMap.h
	egrep -v "#include <pixinth>|#include \"pixintdefs.h\"" < tmp.pixintMap.h > pixintMap.h
	/bin/rm tmp.pixintMap.h

intpixMap.cc intpixMap.h:
	$(GENCLASS) -2 int val Pix val Map intpix
	$(PREPEND_HEADER) intpixMap.h "#include <Pix.h>"
	$(PREPEND_HEADER) intpixMap.h "#include <util/container/pixintdefs.h>"
	/bin/mv intpixMap.h tmp.intpixMap.h
	egrep -v "#include <intpixh>|#include \"intpixdefs.h\"" < tmp.intpixMap.h > intpixMap.h
	/bin/rm tmp.intpixMap.h

pixvpMap.cc pixvpMap.h:
	$(GENCLASS) -2 Pix val VoidPtr val Map pixvp
	$(PREPEND_HEADER) pixvpMap.h "#include <Pix.h>"
	$(PREPEND_HEADER) pixvpMap.h "#include <util/container/voidptr.h>"
	$(PREPEND_HEADER) pixvpMap.h "#include <util/container/pixvpdefs.h>"
	/bin/mv pixvpMap.h tmp.pixvpMap.h
	egrep -v "#include <pixvph>|#include \"pixvpdefs.h\"" < tmp.pixvpMap.h > pixvpMap.h
	/bin/rm tmp.pixvpMap.h

pixvp$(MAPIMPL)Map.cc pixvp$(MAPIMPL)Map.h:
	$(GENCLASS) -2 Pix val VoidPtr val $(MAPIMPL)Map pixvp
	$(PREPEND_HEADER) pixvp$(MAPIMPL)Map.h "#include <Pix.h>"
	$(PREPEND_HEADER) pixvp$(MAPIMPL)Map.h "#include <util/container/voidptr.h>"

pixint$(MAPIMPL)Map.cc pixint$(MAPIMPL)Map.h:
	$(GENCLASS) -2 Pix val int val $(MAPIMPL)Map pixint
	$(PREPEND_HEADER) pixint$(MAPIMPL)Map.h "#include <Pix.h>"

intpix$(MAPIMPL)Map.cc intpix$(MAPIMPL)Map.h:
	$(GENCLASS) -2 int val Pix val $(MAPIMPL)Map intpix
	$(PREPEND_HEADER) intpix$(MAPIMPL)Map.h "#include <Pix.h>"

voidptrSet.cc voidptrSet.h:
	$(GENCLASS) VoidPtr val Set voidptr
	$(PREPEND_HEADER) voidptrSet.h "#include <util/container/voidptr.h>"
	mv voidptrSet.h foo.h
	cat foo.h | sed 's/\"voidptrdefs.h\"/<util\/container\/voidptrdefs.h>/g' > voidptrSet.h
	/bin/rm foo.h

voidptr$(SETIMPL)Set.cc voidptr$(SETIMPL)Set.h:
	$(GENCLASS) VoidPtr val $(SETIMPL)Set voidptr
	$(PREPEND_HEADER) voidptr$(SETIMPL)Set.h "#include <util/container/voidptr.h>"
