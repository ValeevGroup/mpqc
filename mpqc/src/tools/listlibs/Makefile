TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
INCLUDE += -I. -I$(SRCDIR)

#BIN_OR_LIB = BIN
TARGET_TO_MAKE = listlibs

DISTFILES = listlibs.cc Makefile

DEPENDCXXSRC = listlibs.cc

BINOBJ = listlibs.o

LD = $(CXX)

listlibs: $(BINOBJ) $(LIBS)
	$(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS)

install:: listlibs
	$(INSTALL) $(INSTALLBINOPTS) $< $(BINDIR)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(BINOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.o=.d)
endif

