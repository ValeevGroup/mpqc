TOPDIR=../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = mkclasses
BIN_OR_LIB = BIN
LD = $(CXX)

DISTFILES = $(CXXFILES)

CXXFILES = mkclasses.cc

BINOBJ = $(CXXFILES:%.cc=%.o)

DEPENDINCLUDE =

default:: mkclasses

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(BINOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.o=.d)
endif

