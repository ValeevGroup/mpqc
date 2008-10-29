TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
include $(TOPDIR)/lib/Makedirlist
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

ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_CCA),yes)
  DEFINES += -DHAVE_CHEMISTRY_CCA
endif

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

molrender: $(BINOBJ) $(LIBS)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

install:: molrender
	$(MKINSTALLDIRS) $(installroot)$(bindir)
	$(LTINST) $(INSTALL) $(INSTALLBINOPT) $< $(installroot)$(bindir)
	cd $(TOPDIR); CONFIG_FILES=src/bin/molrender/tkmolrender CONFIG_HEADERS= ./config.status
	$(INSTALL) $(INSTALLBINOPT) tkmolrender $(installroot)$(bindir)
	/bin/rm -f tkmolrender

clean::
	/bin/rm -f molrender

distclean::
	/bin/rm -f molrender

$(BINOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.o=.d)
endif

