TOPDIR = ../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libintv2
BIN_OR_LIB = LIB

SGENFLAGS = -ifndef _chemistry_qc_intv2_libint_h

# MG is the maximum angular momentum for which we will use
# the generated build routines.
MG=none
ifeq ($(ARCH),NCUBE_V2)
  MG=2
endif
ifeq ($(ARCH),NCUBE_V3)
  MG=2
endif
ifeq ($(ARCH),SABER)
  MG=0
endif
ifeq ($(MG),none)
  MG=3
endif

ifeq ($(ARCH),NCUBE_V2)
  COPTIONS = -O2
endif
ifeq ($(ARCH),I860)
  COPTIONS = -O3 -Knoieee -Mnoframe
endif

ifeq ($(ARCH),SGI)
  CC = $(GCC)
  COPTIONS = $(GCOPTIONS)
  ifneq ($(CC),$(GCC))
    # this was needed for some SGI's
    #CFLAGS := $(CFLAGS) -Olimit 2000
    CFLAGS := $(CFLAGS_DEBUG)
  endif
endif

CSRC = bounds.c storage.c  buildgc.c  comp_1e.c  comp_0e.c  normalize.c  \
       offsets.c int_fjt.c  init2e.c  atominfo.c  basis.c  comp_erep.c  \
       shiftgc.c int_print.c utils.c

XCSRC = shellprnt.c  shellrbc0.c  shellrcv0.c  shellsbc0.c  shellsnd0.c  \
        shellasgn.c  shellbrd.c  shellbwr.c  shellfrd.c  shellfwr.c \
        shellzero.c

SGENFILE = atoms.sgen

SGENSRC = \
  atomsallc.c atomsbc0.c atomsbrd.c atomsbwr.c atomsfree.c \
  atomsinit.c atomsip.c atomsprnt.c atomsrbc0.c atomssbc0.c \
  atomsfrd.c atomsfwr.c atomsiseq.c atomsasgn.c atomssnd0.c \
  atomsrcv0.c atomszero.c

INC = bounds.h buildinter.h fjttable.h genbuild.h \
               int_flags.h int_libv2.h int_macros.h int_types.h inter.h
SGENINC = $(SGENSRC:%.c=%.h)
GENINC = MG.h $(CSRC:%.c=%.gbl) $(CSRC:%.c=%.lcl)

LIBSRC = $(CSRC) $(XCSRC) $(SGENSRC)
LIBOBJ = $(LIBSRC:%.c=%.$(OBJSUF))

DISTFILES = $(CSRC) $(XCSRC) $(INC) atoms.sgen Makefile LIBS.h

DEPENDINCLUDE = $(INC) $(GENINC) $(SGENINC)

#############################################################

default:: $(DEPENDINCLUDE)
default:: libintv2.a

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

#############################################################

MG.h:
	echo '#define MG $(MG)' > MG.h
