TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
include $(TOPDIR)/lib/Makedirlist
INCLUDE += -I. -I$(SRCDIR)
CXXINCLUDE += -I. -I$(SRCDIR)

TARGET_TO_MAKE = mpqc

CXXSRC = mpqc.cc
CSRC =
INC =

DEPENDINCLUDE = $(INC)

BINOBJ = $(CXXSRC:%.cc=%.o) $(CSRC:%.c=%.o)

DISTFILES = $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h

default:: mpqc

XDEF =
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_PSI),yes)
XDEF += -DHAVE_CHEMISTRY_QC_PSI
endif
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_CC),yes)
XDEF += -DHAVE_CHEMISTRY_QC_CC
endif

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(XDEF) $(SRCDIR)/LIBS.h)

ifeq ($(ARCH),I860)
  DEFINES += -DLINPACKD=-lkmath
endif

#################################################################

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

mpqc: $(BINOBJ) $(LIBS)
	$(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS)

install:: mpqc
	$(INSTALL) $(INSTALLDIROPT) $(bindir)
	$(INSTALL) $(INSTALLBINOPT) $< $(bindir)

distclean::
	/bin/rm -f mpqc

#################################################################

$(BINOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.o=.d)
endif

