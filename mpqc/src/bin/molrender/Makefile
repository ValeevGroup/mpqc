TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
INCLUDE += -I. -I$(SRCDIR)
CXXINCLUDE += -I. -I$(SRCDIR)

TARGET_TO_MAKE = molrender

CXXSRC = main.cc
CSRC =
INC =

DEPENDINCLUDE = $(INC)

BINOBJ = $(CXXSRC:%.cc=%.o) $(CSRC:%.c=%.o)

DISTFILES = $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h molrender.in

default:: molrender

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

molrender: $(BINOBJ) $(LIBS)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

install:: molrender
	$(INSTALL) $(INSTALLDIROPT) $(bindir)
	$(LTINST) $(INSTALL) $(INSTALLBINOPT) $< $(bindir)
	cd $(TOPDIR); CONFIG_FILES=src/bin/molrender/tkmolrender CONFIG_HEADERS= ./config.status
	$(INSTALL) $(INSTALLBINOPT) tkmolrender $(bindir)
	/bin/rm -f tkmolrender

$(BINOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.o=.d)
endif

