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
STACKIMPL = SL

CXXSRC = class.cc proxy.cc

GENCXXSRC = classMap.cc class$(MAPIMPL)Map.cc classkeySet.cc \
            classkey$(SETIMPL)Set.cc
# classStack.cc
#            class$(STACKIMPL)Stack.cc class$(STACKIMPL)List.cc \
#            classl$(STACKIMPL)List.cc

LIBOBJ = $(CXXSRC:%.cc=%.o) $(GENCXXSRC:%.cc=%.o)
GENSRC = $(GENCXXSRC)

CLASSES = DescribedClass DescribedClassProxy

INC = class.h proxy.h
GENINC = classMap.h classImplMap.h classkeySet.h classkeyImplSet.h \
         classRAVLMap.h classkeyAVLSet.h \
         clasmacr.h
# classSLStack.h classStack.h classlSLList.h classSLList.h


DEPENDINCLUDE = $(INC) $(GENINC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libclass

TESTCXXSRC = classtest.cc
DISTFILES = $(CXXSRC) $(INC) Makefile LIBS.h $(TESTCXXSRC)

default:: $(DEPENDINCLUDE)
default:: dlinfo

interface:: $(DEPENDINCLUDE)

clasmacr.h: clastmpl.h $(TMPL2MACDEP)
	$(TMPL2MAC) $< > clasmacr.h

# a test program
LD = $(CXX)
classtest: classtest.o libclass.$(LIBSUF) libcontainer.$(LIBSUF) libref.$(LIBSUF)
	$(LD) $(LDFLAGS) -o classtest $^ $(SYSLIBS)

bug: bug.cc
	$(CXX) $(CXXFLAGS) -o bug $^

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(TESTCXXSRC:%.cc=%.d) $(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(TESTCXXSRC:%.cc=%.d) $(LIBOBJ:.o=.d)
endif

####### Rules for generating source code ######################################

classMap.cc classMap.h:
	$(GENCLASS) -2 ClassKey val ClassDescP val Map class
	$(PREPEND_HEADER) classMap.h "#include <util/class/class.h>"
	$(PREPEND_HEADER) classMap.h "#include <util/class/classdefs.h>"
	/bin/mv classMap.h tmp.classMap.h
	grep -v "#include \"classdefs.h\"" < tmp.classMap.h > classMap.h
	/bin/rm tmp.classMap.h

class$(MAPIMPL)Map.cc class$(MAPIMPL)Map.h:
	$(GENCLASS) -2 ClassKey val ClassDescP val $(MAPIMPL)Map class
	$(PREPEND_HEADER) class$(MAPIMPL)Map.h "#include <util/class/class.h>"

classImplMap.h:
	echo "#include <util/class/class$(MAPIMPL)Map.h>" > classImplMap.h
	echo "#define MAPCTOR ClassKeyClassDescP"$(MAPIMPL)"Map(0)" >> classImplMap.h

classkeyPlex.cc classkeyPlex.h:
	$(GENCLASS) ClassKey val Plex classkey
	$(PREPEND_HEADER) classkeyPlex.h "#include <util/class/class.h>"

classkeyXPlex.cc classkeyXPlex.h:
	$(GENCLASS) ClassKey val XPlex classkey
	$(PREPEND_HEADER) classkeyXPlex.h "#include <util/class/class.h>"

classkeySet.cc classkeySet.h:
	$(GENCLASS) ClassKey val Set classkey
	$(PREPEND_HEADER) classkeySet.h "#include <util/class/class.h>"
	$(PREPEND_HEADER) classkeySet.h "#include <util/class/classkeydefs.h>"
	/bin/mv classkeySet.h tmp.classkeySet.h
	grep -v "#include \"classkeydefs.h\"" < tmp.classkeySet.h > classkeySet.h
	/bin/rm tmp.classkeySet.h

classkey$(SETIMPL)Set.cc classkey$(SETIMPL)Set.h:
	$(GENCLASS) ClassKey val $(SETIMPL)Set classkey
	$(PREPEND_HEADER) classkey$(SETIMPL)Set.h "#include <util/class/class.h>"

classkeyImplSet.h:
	echo "#include <util/class/classkey$(SETIMPL)Set.h>" > classkeyImplSet.h
	echo "#define SETCTOR ClassKey"$(SETIMPL)"Set()" >> classkeyImplSet.h

classStack.cc classStack.h:
	$(GENCLASS) ClassDescP val Stack class
	$(PREPEND_HEADER) classStack.h "#include <util/class/class.h>"

class$(STACKIMPL)Stack.cc class$(STACKIMPL)Stack.h:
	$(GENCLASS) ClassDescP val $(STACKIMPL)Stack class
	$(PREPEND_HEADER) class$(STACKIMPL)Stack.h "#include <util/class/class.h>"

classSLList.cc classSLList.h:
	$(GENCLASS) ClassDescP val SLList class
	$(PREPEND_HEADER) classSLList.h "#include <util/class/class.h>"

classlSLList.cc classlSLList.h:
	$(GENCLASS) ClassDescPSLList val SLList classl
	$(PREPEND_HEADER) classlSLList.h "#include <util/class/classSLList.h>"
	$(PREPEND_HEADER) classlSLList.h "#include <util/class/class.h>"

classImplStack.h:
	echo "#include <util/class/class$(STACKIMPL)Stack.h" > classImplStack.h
	echo "#define STACKCTOR ClassDescP"$(STACKIMPL)"Stack()" >> classImplStack.h

