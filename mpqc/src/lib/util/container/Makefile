TOPDIR=../../..
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

REALCXXSRC = array.cc set.cc ref.cc voidptr.cc gnuavlse.cc

GENSRC = $(GENCXXSRC)

CXXSRC = $(GENCXXSRC) $(REALCXXSRC)
TESTCXX = settest.cc settestt.cc reftest.cc reftestx.cc

LIBOBJ = $(CXXSRC:%.cc=%.o)

INC = pixmap.h pixvpdefs.h ptrset.h voidptr.h voidptrdefs.h \
      set.h array.h ref.h reftestx.h artem.h ar2tem.h ssarmac.h ssar2mac.h \
      gnuset.h gnuavlse.h pixintdefs.h

GENINC = $(GENCXXSRC:%.cc=%.h) setmacr.h asetmacr.h refmacr.h \
         armac.h ar2mac.h ssarmac.h ssar2mac.h

DEPENDINCLUDE = $(INC) $(GENINC)

DISTFILES = Makefile $(INC) LIBS.h $(REALCXXSRC) $(TESTCXX)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libcontainer

default:: $(DEPENDINCLUDE)
default:: $(TOPDIR)/lib/libcontainer.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libcontainer.so
endif

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

refmacr.h: reftmpl.h $(TMPL2MACDEP)
	$(TMPL2MAC) $< > refmacr.h

settest: settest.o libcontainer.$(LIBSUF)
	$(CXX) $(LDFLAGS) -o settest $^

gnutest: gnutest.o libcontainer.$(LIBSUF)
	$(CXX) $(LDFLAGS) -o gnutest $^

artest: artest.o libcontainer.$(LIBSUF)
	$(CXX) $(LDFLAGS) -o artest $^

settestt: settestt.o libcontainer.$(LIBSUF)
	$(CXX) $(LDFLAGS) -o settestt $^

reftest: reftest.o reftestx.o libcontainer.$(LIBSUF) libmisc.$(LIBSUF)
	$(CXX) $(LDFLAGS) -o reftest $^

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

iclean::
	-rm -f foo.h

distclean::
	-/bin/rm -f $(GENINC)

install_inc:: $(INC) $(GENINC)
	$(INCINSTALL) $^ $(INCDIR)/SC/util/container

$(TESTCXX:%.cc=%.d) $(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

pixintMap.cc pixintMap.h:
	$(GENCLASS) -2 Pix val int val Map pixint
	prepend-header pixintMap.h "#include <Pix.h>"
	prepend-header pixintMap.h "#include <util/container/pixintdefs.h>"
	/bin/mv pixintMap.h tmp.pixintMap.h
	egrep -v "#include <pixinth>|#include \"pixintdefs.h\"" < tmp.pixintMap.h > pixintMap.h
	/bin/rm tmp.pixintMap.h

intpixMap.cc intpixMap.h:
	$(GENCLASS) -2 int val Pix val Map intpix
	prepend-header intpixMap.h "#include <Pix.h>"
	prepend-header intpixMap.h "#include <util/container/pixintdefs.h>"
	/bin/mv intpixMap.h tmp.intpixMap.h
	egrep -v "#include <intpixh>|#include \"intpixdefs.h\"" < tmp.intpixMap.h > intpixMap.h
	/bin/rm tmp.intpixMap.h

pixvpMap.cc pixvpMap.h:
	$(GENCLASS) -2 Pix val VoidPtr val Map pixvp
	prepend-header pixvpMap.h "#include <Pix.h>"
	prepend-header pixvpMap.h "#include <util/container/voidptr.h>"
	prepend-header pixvpMap.h "#include <util/container/pixvpdefs.h>"
	/bin/mv pixvpMap.h tmp.pixvpMap.h
	egrep -v "#include <pixvph>|#include \"pixvpdefs.h\"" < tmp.pixvpMap.h > pixvpMap.h
	/bin/rm tmp.pixvpMap.h

pixvp$(MAPIMPL)Map.cc pixvp$(MAPIMPL)Map.h:
	$(GENCLASS) -2 Pix val VoidPtr val $(MAPIMPL)Map pixvp
	prepend-header pixvp$(MAPIMPL)Map.h "#include <Pix.h>"
	prepend-header pixvp$(MAPIMPL)Map.h "#include <util/container/voidptr.h>"

pixint$(MAPIMPL)Map.cc pixint$(MAPIMPL)Map.h:
	$(GENCLASS) -2 Pix val int val $(MAPIMPL)Map pixint
	prepend-header pixint$(MAPIMPL)Map.h "#include <Pix.h>"

intpix$(MAPIMPL)Map.cc intpix$(MAPIMPL)Map.h:
	$(GENCLASS) -2 int val Pix val $(MAPIMPL)Map intpix
	prepend-header intpix$(MAPIMPL)Map.h "#include <Pix.h>"

voidptrSet.cc voidptrSet.h:
	$(GENCLASS) VoidPtr val Set voidptr
	prepend-header voidptrSet.h "#include <util/container/voidptr.h>"
	mv voidptrSet.h foo.h
	cat foo.h | sed 's/\"voidptrdefs.h\"/<util\/container\/voidptrdefs.h>/g' > voidptrSet.h
	/bin/rm foo.h

voidptr$(SETIMPL)Set.cc voidptr$(SETIMPL)Set.h:
	$(GENCLASS) VoidPtr val $(SETIMPL)Set voidptr
	prepend-header voidptr$(SETIMPL)Set.h "#include <util/container/voidptr.h>"
