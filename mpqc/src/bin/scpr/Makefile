TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
include $(TOPDIR)/lib/Makedirlist
INCLUDE += -I. -I$(SRCDIR)
CXXINCLUDE += -I. -I$(SRCDIR)

TARGET_TO_MAKE = scpr

CXXSRC = scpr.cc
CSRC =
INC = 

DEPENDINCLUDE = $(INC)

BINOBJ = $(CXXSRC:%.cc=%.$(OBJSUF)) $(CSRC:%.c=%.$(OBJSUF))

DISTFILES = $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h

default:: scpr

XDEF=
XLIBS=
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_CINTS),yes)
 ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_MBPTR12),yes)
 XDEF += -DHAVE_CHEMISTRY_QC_MBPTR12
 endif
XDEF += -DHAVE_CHEMISTRY_QC_CINTS
endif

ifeq ($(HAVE_LIBINT),yes)
  XLIBS += -lint
endif
ifeq ($(HAVE_LIBR12),yes)
  XLIBS += -lr12
endif
ifeq ($(HAVE_LIBDERIV),yes)
  XLIBS += -lderiv
endif

ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_CCA),yes)
  XDEF += -DHAVE_CHEMISTRY_CCA
endif

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(XDEF) $(SRCDIR)/LIBS.h)

#################################################################

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

scpr: $(BINOBJ) $(LIBS)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(XLIBS) $(SYSLIBS) $(LTLINKBINOPTS)

install:: scpr
	$(MKINSTALLDIRS) $(installroot)$(bindir)
	$(LTINST) $(INSTALL) $(INSTALLBINOPT) $< $(installroot)$(bindir)

clean::
	/bin/rm -f scpr

distclean::
	/bin/rm -f scpr

#################################################################

$(BINOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.$(OBJSUF)=.d)
endif
