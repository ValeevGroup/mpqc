TOPDIR=../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

LD = $(CXX)

TARGET_TO_MAKE = mpqcic
#BIN_OR_LIB = BIN

CXXSRC = mpqcic.cc mpqc_int.cc esp.cc mp2.cc
CSRC = gcollect.c
INC = mpqc_int.h

DEPENDINCLUDE = $(INC)

BINOBJ = $(CXXSRC:%.cc=%.o) $(CSRC:%.c=%.o)

DISTFILES = $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h

default:: mpqcic

MPQCLIBSDEFINE=
DEFAULT = mpqcicproc

ifeq ($(ARCH),SGI)
  DEFAULT = mpqcicproc
endif

ifeq ($(ARCH),PARAGON)
  MPQCLIBSDEFINE=-DUSE_PICLPARAGON
  DEFAULT=mpqcic
endif

ifeq ($(ARCH),I860)
  MPQCLIBSDEFINE=-DUSE_PICLNODE
  DEFAULT=mpqcic
endif

ifeq ($(ARCH),NCUBE_V2)
  MPQCLIBSDEFINE=-DUSE_PICLN6400
  DEFAULT = mpqcic
endif

ifeq ($(ARCH),NCUBE_V3)
  MPQCLIBSDEFINE=-DUSE_PICLN6400
  DEFAULT = mpqcic
endif


LIBS := $(shell $(LISTLIBS) $(MPQCLIBSDEFINE) $(INCLUDE) $(SRCDIR)/LIBS.h)
#PROCLIBS = $(shell $(LISTLIBS) -DUSE_PICLPROC $(INCLUDE) $(SRCDIR)/LIBS.h)
#PVMLIBS = $(shell $(LISTLIBS) -DUSE_PICLPVM $(INCLUDE) $(SRCDIR)/LIBS.h)
#SHMLIBS = $(shell $(LISTLIBS) -DUSE_PICLSHM $(INCLUDE) $(SRCDIR)/LIBS.h)


ifeq ($(ARCH),NCUBE_V2)
  CFLAGS := $(CFLAGS) -Ncomm  500000 -Nheap 2500000 -Nstack 64000 -d 1
  DATE := -DCOMPILE_DATE=NULL
endif

ifeq ($(ARCH),NCUBE_V3)
  CFLAGS := $(CFLAGS) -Ncomm  500000 -Nheap 2500000 -Nstack 64000 -d 1
endif

ifeq ($(ARCH),I860)
  LIBS := $(shell $(LISTLIBS) $(MPQCLIBSDEFINE) $(INCLUDE) $(SRCDIR)/LIBS.h)
  DEFINES += -DLINPACKD=-lkmath
endif

ifeq ($(ARCH),SGI)
  CC = $(GCC)
  COPTIONS = $(GCOPTIONS_DEBUG)
endif

#################################################################

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules


mpqcicproc: $(BINOBJ) $(LIBS)
	$(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS)

install_mpqcicproc:: mpqcicproc
	$(INSTALL) $(INSTALLBINOPTS) $< $(BINDIR)

install_mpqcicpvm:: mpqcicpvm
	$(INSTALL) $(INSTALLBINOPTS) $< $(BINDIR)
mpqcicpvm: $(BINOBJ) \
	  $(shell $(LISTLIBS) -DUSE_PICLPVM $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS)

install_mpqcicshm:: mpqcicshm
	$(INSTALL) $(INSTALLBINOPTS) $< $(BINDIR)
mpqcicshm: $(BINOBJ) \
	  $(shell $(LISTLIBS) -DUSE_PICLSHM $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS)

install_mpqcic:: mpqcic
	$(INSTALL) $(INSTALLBINOPTS) $< $(BINDIR)

mpqcic: $(BINOBJ) $(LIBS)
	$(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS)


#################################################################

$(BINOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.o=.d)
endif

