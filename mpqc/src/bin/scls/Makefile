TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
INCLUDE += -I. -I$(SRCDIR)
CXXINCLUDE += -I. -I$(SRCDIR)

TARGET_TO_MAKE = scls

CXXSRC = scls.cc
CSRC =
INC = 

DEPENDINCLUDE = $(INC)

BINOBJ = $(CXXSRC:%.cc=%.o) $(CSRC:%.c=%.o)

DISTFILES = $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h

default:: scls

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

#################################################################

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

scls: $(BINOBJ) $(LIBS)
	$(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS)

install:: scls
	$(INSTALL) $(INSTALLDIROPT) $(bindir)
	$(INSTALL) $(INSTALLBINOPT) $< $(bindir)

distclean::
	/bin/rm -f scls

#################################################################

$(BINOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.o=.d)
endif
