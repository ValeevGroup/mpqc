TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

CXXSRC = symop.cc rep.cc irrep.cc chartab.cc pointgrp.cc \
         maketab.cc tetra.cc ico.cc

LIBOBJ = $(CXXSRC:%.cc=%.o)

INC =  pointgrp.h

DEPENDINCLUDE = $(INC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libsymmetry

DISTFILES = $(CXXSRC) $(INC) Makefile LIBS.h

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
testpg: testpg.o $(LIBS)
	$(CXX) $^ -o $@ $(SYSLIBS)

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

