TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

REALCXXSRC = ref.cc identity.cc

GENCXXSRC =
GENSRC = $(GENCXXSRC)

CXXSRC = $(GENCXXSRC) $(REALCXXSRC)
TESTCXX = reftest.cc reftestx.cc

LIBOBJ = $(CXXSRC:%.cc=%.o)

INC = ref.h reftestx.h reftmpl.h identity.h

GENINC = refmacr.h

DEPENDINCLUDE = $(INC) $(GENINC)

DISTFILES = Makefile $(INC) LIBS.h $(REALCXXSRC) $(TESTCXX)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libref

default:: $(DEPENDINCLUDE)
default:: $(TOPDIR)/lib/libref.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libref.so
endif

interface:: $(DEPENDINCLUDE)

refmacr.h: reftmpl.h $(TMPL2MACDEP)
	$(TMPL2MAC) $< > refmacr.h

reftest: reftest.o reftestx.o libref.$(LIBSUF)
	$(CXX) $(LDFLAGS) -o reftest $^

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

iclean::
	-rm -f foo.h

distclean::
	-/bin/rm -f $(GENINC)

$(TESTCXX:%.cc=%.d) $(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif
