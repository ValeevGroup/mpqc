TOPDIR = ../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libscf
BIN_OR_LIB = LIB

SRCS = accum.cc grscf.cc effh.cc hcore.cc

CLASSES = GRSCF

LIBOBJ= $(SRCS:%.cc=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)
default:: dlinfo
default:: $(TOPDIR)/lib/libscf.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libscf.so
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

