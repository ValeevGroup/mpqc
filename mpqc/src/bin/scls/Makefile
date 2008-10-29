TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
include $(TOPDIR)/lib/Makedirlist
INCLUDE += -I. -I$(SRCDIR)
CXXINCLUDE += -I. -I$(SRCDIR)

TARGET_TO_MAKE = scls

CXXSRC = scls.cc
CSRC =
INC = 

DEPENDINCLUDE = $(INC)

BINOBJ = $(CXXSRC:%.cc=%.$(OBJSUF)) $(CSRC:%.c=%.$(OBJSUF))

DISTFILES = $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h

default:: scls

ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_CCA),yes)
  XDEF += -DHAVE_CHEMISTRY_CCA
endif

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(XDEF) $(SRCDIR)/LIBS.h)

#################################################################

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

scls: $(BINOBJ) $(LIBS)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

install:: scls
	$(MKINSTALLDIRS) $(installroot)$(bindir)
	$(LTINST) $(INSTALL) $(INSTALLBINOPT) $< $(installroot)$(bindir)

clean::
	/bin/rm -f scls

distclean::
	/bin/rm -f scls

#################################################################

$(BINOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.$(OBJSUF)=.d)
endif
