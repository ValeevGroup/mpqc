TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

#BIN_OR_LIB = BIN
TARGET_TO_MAKE = mpqc++
#LD = $(CXX)

DISTFILES = Makefile main.cc

BINOBJ = main.o

default:: mpqc++

install:: mpqc++
	$(INSTALL) $(INSTALLBINOPTS) $< $(BINDIR)

mpqc++: $(BINOBJ) \
	 $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SYSLIBS)


include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(BINOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.o=.d)
endif

