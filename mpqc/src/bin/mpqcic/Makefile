TOPDIR=../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
INCLUDE += -I. -I$(SRCDIR)
CXXINCLUDE += -I. -I$(SRCDIR)

TARGET_TO_MAKE = mpqcic

CXXSRC = mpqcic.cc mpqc_int.cc
CSRC =
INC = mpqc_int.h g92_int.h

DEPENDINCLUDE = $(INC)

BINOBJ = $(CXXSRC:%.cc=%.o) $(CSRC:%.c=%.o)

DISTFILES = $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h

default:: mpqcic

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

ifeq ($(ARCH),NCUBE_V2)
  CFLAGS := $(CFLAGS) -Ncomm  500000 -Nheap 2500000 -Nstack 64000 -d 1
  DATE := -DCOMPILE_DATE=NULL
endif

ifeq ($(ARCH),NCUBE_V3)
  CFLAGS := $(CFLAGS) -Ncomm  500000 -Nheap 2500000 -Nstack 64000 -d 1
endif

ifeq ($(ARCH),I860)
  DEFINES += -DLINPACKD=-lkmath
endif

#################################################################

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

install_mpqcic:: mpqcic
	$(INSTALL) $(INSTALLBINOPTS) $< $(BINDIR)

mpqcic: $(BINOBJ) $(LIBS)
	$(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS)

mpqcg92: g92.o run_g92.o g92freq.o mpqc_int.o $(LIBS)
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

