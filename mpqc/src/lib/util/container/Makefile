TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

SETIMPL = AVL
MAPIMPL = RAVL

GENCXXSRC = pixvpMap.cc pixvp$(MAPIMPL)Map.cc \
            voidptrSet.cc voidptr$(SETIMPL)Set.cc

REALCXXSRC = array.cc set.cc ref.cc

GENSRC = $(GENCXXSRC)

CXXSRC = $(GENCXXSRC) $(REALCXXSRC)

LIBOBJ = $(CXXSRC:%.cc=%.o)

INC = pixmap.h pixvpdefs.h ptrset.h voidptr.h voidptrdefs.h \
      set.h array.h

GENINC = $(GENCXXSRC:%.cc=%.h)

DEPENDINCLUDE = $(INC) $(GENINC)

DISTFILES = Makefile $(INC) LIBS.h settest.h $(REALCXXSRC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libcontainer

default:: $(DEPENDINCLUDE)
default:: libcontainer.a

settest: settest.o libcontainer.a
	$(CXX) $(CXXFLAGS) -o settest $^

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

iclean::
	-rm -f foo.h

install_inc:: $(INC) $(GENINC)
	$(INCINSTALL) $^ $(INCDIR)/SC/util/container

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

pixvpMap.cc pixvpMap.h:
	genclass -2 Pix val VoidPtr val Map pixvp
	prepend-header pixvpMap.h "#include <Pix.h>"
	prepend-header pixvpMap.h "#include <util/container/voidptr.h>"
	prepend-header pixvpMap.h "#include <util/container/pixvpdefs.h>"
	/bin/mv pixvpMap.h tmp.pixvpMap.h
	egrep -v "#include <pixvph>|#include \"pixvpdefs.h\"" < tmp.pixvpMap.h > pixvpMap.h
	/bin/rm tmp.pixvpMap.h

pixvp$(MAPIMPL)Map.cc pixvp$(MAPIMPL)Map.h:
	genclass -2 Pix val VoidPtr val $(MAPIMPL)Map pixvp
	prepend-header pixvp$(MAPIMPL)Map.h "#include <Pix.h>"
	prepend-header pixvp$(MAPIMPL)Map.h "#include <util/container/voidptr.h>"

voidptrSet.cc voidptrSet.h:
	genclass VoidPtr val Set voidptr
	prepend-header voidptrSet.h "#include <util/container/voidptr.h>"
	mv voidptrSet.h foo.h
	cat foo.h | sed 's/\"voidptrdefs.h\"/<util\/container\/voidptrdefs.h>/g' > voidptrSet.h
	/bin/rm foo.h

voidptr$(SETIMPL)Set.cc voidptr$(SETIMPL)Set.h:
	genclass VoidPtr val $(SETIMPL)Set voidptr
	prepend-header voidptr$(SETIMPL)Set.h "#include <util/container/voidptr.h>"
