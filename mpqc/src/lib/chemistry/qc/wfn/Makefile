TOPDIR = ../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

#CXXFLAGS = $(CXXFLAGS_DEBUG)

TARGET_TO_MAKE = libwfn
BIN_OR_LIB = LIB

SRCS = wfn.cc densval.cc natorbval.cc obwfn.cc density.cc hcorewfn.cc

CLASSES = ElectronDensity OneBodyWavefunction Wavefunction

LIBOBJ= $(SRCS:%.cc=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)
default:: dlinfo
default:: $(TOPDIR)/lib/libwfn.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libwfn.so
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

