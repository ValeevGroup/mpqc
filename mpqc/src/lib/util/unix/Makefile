TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

CSRC = cutils.c

CXXSRC =

ifeq ($(ARCH),none)
  CXXSRC += daemon.cc usignal.cc
endif

LIBOBJ = $(CSRC:%.c=%.o) $(CXXSRC:%.cc=%.o)

INC =  cct_cprot.h daemon.h usignal.h

DEPENDINCLUDE = $(INC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libunix

DISTFILES = $(CSRC) $(CXXSRC) $(INC) Makefile LIBS.h

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

