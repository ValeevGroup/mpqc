TOPDIR=../../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libintv3
BIN_OR_LIB = LIB

SGENFLAGS = -ifndef _chemistry_qc_intv3_libint_h

TESTCSRC = inttest.cc

CXXSRC = fjt.cc offsets.cc tformv3.cc \
         int1e.cc comp1e.cc \
         int2e.cc comp2e.cc build2e.cc print2e.cc shift2e.cc init2e.cc \
         comp2e3c.cc \
         bounds.cc store.cc storage.cc \
         obintv3.cc tbintv3.cc intv3.cc \
         array.cc

INC = fjt.h int1e.h int2e.h utils.h tformv3.h \
      tbintv3.h obintv3.h intv3.h cartitv3.h

LIBSRC = $(CXXSRC) $(CSRC) $(XCSRC) $(SGENSRC)
LIBOBJ = $(LIBSRC:%.c=%.$(OBJSUF))
LIBOBJ := $(LIBOBJ:%.cc=%.$(OBJSUF))

DISTFILES = $(TESTCSRC) $(INC) atoms.sgen Makefile LIBS.h

DEPENDINCLUDE = $(INC) $(GENINC) $(SGENINC)

#############################################################

default:: $(DEPENDINCLUDE)
default:: $(TOPDIR)/lib/libintv3.a
interface:: $(DEPENDINCLUDE)
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libintv3.so
endif

interface:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

inttest: inttest.o \
	$(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(CXX) $(LDFLAGS) -o inttest $^ $(SYSLIBS)

inttest.o: inttest.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

$(TESTOBJ:.o=.d) $(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTOBJ:.o=.d)
endif

#############################################################
