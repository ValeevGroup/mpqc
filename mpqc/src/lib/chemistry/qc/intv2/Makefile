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

TESTCSRC = inttest.cc

CSRC = bounds.c buildgc.c  comp_1e.c  comp_0e.c  normalize.c  \
       offsets.c int_fjt.c  init2e.c  atominfo.c  basis.c  comp_erep.c  \
       shiftgc.c int_print.c utils.c comp_erep23.c transform.c

CXXSRC = basisk.cc storage.cc

XCSRC = shellprnt.c  shellrbc0.c  shellrcv0.c  shellsbc0.c  shellsnd0.c  \
        shellasgn.c  shellzero.c

SGENFILE = atoms.sgen

SGENSRC = \
  atomsallc.c atomsbc0.c atomsfree.c \
  atomsinit.c atomsprnt.c atomsrbc0.c atomssbc0.c \
  atomsiseq.c atomsasgn.c atomssnd0.c \
  atomsrcv0.c atomszero.c

INC = bounds.h buildinter.h fjttable.h genbuild.h transform.h \
      storage.h int_flags.h int_libv2.h int_macros.h int_types.h inter.h
SGENINC = $(SGENSRC:%.c=%.h)
GENINC = MG.h $(CSRC:%.c=%.gbl) $(CSRC:%.c=%.lcl)

LIBSRC = $(CSRC) $(XCSRC) $(SGENSRC) $(CXXSRC)
LIBOBJ = $(LIBSRC:%.c=%.$(OBJSUF))
LIBOBJ := $(LIBOBJ:%.cc=%.$(OBJSUF))

DISTFILES = $(TESTCSRC) $(CSRC) $(XCSRC) $(INC) atoms.sgen Makefile LIBS.h

DEPENDINCLUDE = $(INC) $(GENINC) $(SGENINC)

#############################################################

default:: $(DEPENDINCLUDE)
default:: $(TOPDIR)/lib/libintv2.a
interface:: $(DEPENDINCLUDE)
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libintv2.so
endif

interface:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

inttest: inttest.o \
	$(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(CXX) $(LDFLAGS) -o inttest $^ $(SYSLIBS)

inttest.o: inttest.cc
	$(CXX) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

$(TESTOBJ:.o=.d) $(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTOBJ:.o=.d)
endif

#############################################################

MG.h:
	echo '#define MG $(MG)' > MG.h

#buildgc.o: buildgc.c
#	$(CC) $(CFLAGS) $(HIGHOPT_CFLAGS) $(NOALIAS_CFLAGS) -c $<
