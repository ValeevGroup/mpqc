TOPDIR=../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

LD = $(CXX)

TARGET_TO_MAKE = mpqcic

CXXSRC = mpqcic.cc mpqc_int.cc opt2_v1.cc opt2_v2.cc
CSRC = scf_ffo.c opt2_fock.c gcollect.c
INC = mpqc_int.h opt2.h g92_int.h

DEPENDINCLUDE = $(INC) \
                scf_ffo.gbl scf_ffo.lcl opt2_fock.gbl opt2_fock.lcl

BINOBJ = $(CXXSRC:%.cc=%.o) $(CSRC:%.c=%.o)

DISTFILES = $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h

default:: mpqcic

MPQCLIBSDEFINE=
DEFAULT = mpqcicproc

ifeq ($(ARCH),SGI)
  DEFAULT = mpqcicproc
endif

ifeq ($(ARCH),SUNMOS)
  MPQCLIBSDEFINE=-DUSE_PICLPARAGON
  DEFAULT=mpqcic
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

g92: g92.o run_g92.o g92freq.o mpqc_int.o $(LIBS)
	$(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS)

writeavs: writeavs.o $(LIBS)
	$(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS)

calcesp: calcesp.o $(LIBS)
	$(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS)


#################################################################

$(BINOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.o=.d)
endif

