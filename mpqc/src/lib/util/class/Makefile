TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

SETIMPL = AVL
MAPIMPL = RAVL
STACKIMPL = SL

CXXSRC = class.cc

GENCXXSRC = classMap.cc class$(MAPIMPL)Map.cc classkeySet.cc \
            classkey$(SETIMPL)Set.cc
# classStack.cc
#            class$(STACKIMPL)Stack.cc class$(STACKIMPL)List.cc \
#            classl$(STACKIMPL)List.cc

LIBOBJ = $(CXXSRC:%.cc=%.o) $(GENCXXSRC:%.cc=%.o)
GENSRC = $(GENCXXSRC)

CLASSES = DescribedClass

INC = class.h
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
default:: $(TOPDIR)/lib/libclass.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libclass.so
endif

interface:: $(DEPENDINCLUDE)

clasmacr.h: clastmpl.h $(TMPL2MACDEP)
	$(TMPL2MAC) $< > clasmacr.h

# a test program
LD = $(CXX)
classtest: classtest.o libclass.$(LIBSUF) libcontainer.$(LIBSUF) libmisc.$(LIBSUF)
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
	prepend-header classMap.h "#include <util/class/class.h>"
	prepend-header classMap.h "#include <util/class/classdefs.h>"
	/bin/mv classMap.h tmp.classMap.h
	grep -v "#include \"classdefs.h\"" < tmp.classMap.h > classMap.h
	/bin/rm tmp.classMap.h

class$(MAPIMPL)Map.cc class$(MAPIMPL)Map.h:
	$(GENCLASS) -2 ClassKey val ClassDescP val $(MAPIMPL)Map class
	prepend-header class$(MAPIMPL)Map.h "#include <util/class/class.h>"

classImplMap.h:
	echo "#include "\"class$(MAPIMPL)Map.h\" > classImplMap.h
	echo "#define MAPCTOR ClassKeyClassDescP"$(MAPIMPL)"Map(0)" >> classImplMap.h

classkeyPlex.cc classkeyPlex.h:
	$(GENCLASS) ClassKey val Plex classkey
	prepend-header classkeyPlex.h "#include \"class.h\""

classkeyXPlex.cc classkeyXPlex.h:
	$(GENCLASS) ClassKey val XPlex classkey
	prepend-header classkeyXPlex.h "#include \"class.h\""

classkeySet.cc classkeySet.h:
	$(GENCLASS) ClassKey val Set classkey
	prepend-header classkeySet.h "#include \"class.h\""

classkey$(SETIMPL)Set.cc classkey$(SETIMPL)Set.h:
	$(GENCLASS) ClassKey val $(SETIMPL)Set classkey
	prepend-header classkey$(SETIMPL)Set.h "#include \"class.h\""

classkeyImplSet.h:
	echo "#include "\"classkey$(SETIMPL)Set.h\" > classkeyImplSet.h
	echo "#define SETCTOR ClassKey"$(SETIMPL)"Set()" >> classkeyImplSet.h

classStack.cc classStack.h:
	$(GENCLASS) ClassDescP val Stack class
	prepend-header classStack.h "#include \"class.h\""

class$(STACKIMPL)Stack.cc class$(STACKIMPL)Stack.h:
	$(GENCLASS) ClassDescP val $(STACKIMPL)Stack class
	prepend-header class$(STACKIMPL)Stack.h "#include \"class.h\""

classSLList.cc classSLList.h:
	$(GENCLASS) ClassDescP val SLList class
	prepend-header classSLList.h "#include \"class.h\""

classlSLList.cc classlSLList.h:
	$(GENCLASS) ClassDescPSLList val SLList classl
	prepend-header classlSLList.h "#include \"classSLList.h\""
	prepend-header classlSLList.h "#include \"class.h\""

classImplStack.h:
	echo "#include "\"class$(STACKIMPL)Stack.h\" > classImplStack.h
	echo "#define STACKCTOR ClassDescP"$(STACKIMPL)"Stack()" >> classImplStack.h

