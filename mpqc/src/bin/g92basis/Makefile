TOPDIR=../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

#BIN_OR_LIB = BIN
TARGET_TO_MAKE = g92basis
#LD = $(CXX)

DISTFILES = Makefile g92basis.cc

BINOBJ = g92basis.o

default:: g92basis

install:: g92basis
	$(INSTALL) $(INSTALLBINOPTS) $< $(BINDIR)

g92basis: $(BINOBJ) \
	 $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SYSLIBS)


include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(BINOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.o=.d)
endif

