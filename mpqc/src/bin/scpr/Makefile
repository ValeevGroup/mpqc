TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
INCLUDE += -I. -I$(SRCDIR)
CXXINCLUDE += -I. -I$(SRCDIR)

TARGET_TO_MAKE = scpr

CXXSRC = scpr.cc
CSRC =
INC = 

DEPENDINCLUDE = $(INC)

BINOBJ = $(CXXSRC:%.cc=%.o) $(CSRC:%.c=%.o)

DISTFILES = $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h

default:: scpr

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

#################################################################

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

scpr: $(BINOBJ) $(LIBS)
	$(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS)

install: scpr
	$(INSTALL) $(INSTALLDIROPT) $(bindir)
	$(INSTALL) $(INSTALLBINOPT) $< $(bindir)

distclean::
	/bin/rm -f scpr

#################################################################

$(BINOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.o=.d)
endif
