TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

BIN_OR_LIB = none
TARGET_TO_MAKE = mpqc

#DEFINES += -DUSE_DEBUG

ifeq ($(ARCH),SUN4)
  CFLAGS := $(CFLAGS_DEBUG)
endif

ifeq ($(ARCH),I860)
  PROFLIBDIR = /usr/local/src/parallel/ssd/i860/lib-coff/prof
#  PROFLIB = $(PROFLIBDIR)/libnode.a $(PROFLIBDIR)/libprof.a
#  COPTIONS = -O3 -node -Mnoframe -Knoieee
#  SYSLIBS = -lX11 -lf $(LIBDIR)/libmsqrt.a $(PROFLIB)
#  SYSLIBS = -lX11 $(LIBDIR)/libfnew.a $(LIBDIR)/libmnew.a $(PROFLIB)
  SYSLIBS = -lX11 -lf
endif

MPQCLIBSDEFINE=none
DEFAULT=mpqcproc
ifeq ($(ARCH),SGI)
DEFAULT = mpqcproc
endif
ifeq ($(ARCH),SUN4)
DEFAULT = mpqcproc
endif
ifeq ($(ARCH),PARAGON)
MPQCLIBSDEFINE=-DUSE_PICLPARAGON
DEFAULT = mpqc
endif
ifeq ($(ARCH),NCUBE_V2)
MPQCLIBSDEFINE=-DUSE_PICLN6400
DEFAULT = mpqc
endif
ifeq ($(ARCH),NCUBE_V3)
MPQCLIBSDEFINE=-DUSE_PICLN6400
DEFAULT = mpqc
endif
ifeq ($(ARCH),I860)
MPQCLIBSDEFINE=-DUSE_PICLNODE
DEFAULT = mpqc
endif

default:: $(DEFAULT)

DATE := -DCOMPILE_DATE="\"$(shell date)\""

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

DISTFILES = Makefile mpqcnode.c compdate.c includes.h disclaimer.h

install:: install_$(DEFAULT)


install_mpqcproc:: mpqcproc
	$(INSTALL) $(INSTALLBINOPTS) $< $(BINDIR)
mpqcproc: mpqcnode.$(OBJSUF) compdate.$(OBJSUF) \
	 $(shell $(LISTLIBS) -DUSE_PICLPROC $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(CC) $(CFLAGS) -o $@ $^ $(SYSLIBS)

install_mpqcpvm:: mpqcpvm
	$(INSTALL) $(INSTALLBINOPTS) $< $(BINDIR)
mpqcpvm: mpqcnode.$(OBJSUF) compdate.$(OBJSUF) \
	  $(shell $(LISTLIBS) -DUSE_PICLPVM $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(CC) $(CFLAGS) -o $@ $^ $(SYSLIBS)

install_mpqcshm:: mpqcshm
	$(INSTALL) $(INSTALLBINOPTS) $< $(BINDIR)
mpqcshm: mpqcnode.$(OBJSUF) compdate.$(OBJSUF) \
	  $(shell $(LISTLIBS) -DUSE_PICLSHM $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(CC) $(CFLAGS) -o $@ $^ $(SYSLIBS)

ifneq ($(MPQCLIBSDEFINE),none)
install_mpqc:: mpqc
	$(INSTALL) $(INSTALLBINOPTS) $< $(BINDIR)
mpqc: mpqcnode.$(OBJSUF) compdate.$(OBJSUF) \
	 $(shell $(LISTLIBS) $(MPQCLIBSDEFINE) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS)
endif

compdate.o: compdate.c always_needs_made
	$(CC) $(DEFINES) $(DATE) -c -o compdate.$(OBJSUF) $<
always_needs_made:

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

values::
	@echo PROCLIBS=$(shell $(LISTLIBS) -DUSE_PICLPROC $(INCLUDE) $(SRCDIR)/LIBS.h)

OBJ = mpqcnode.o
$(OBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(OBJ:.o=.d)
endif

