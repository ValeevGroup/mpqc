TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

#CXXFLAGS = $(CXXFLAGS_DEBUG)

CXXSRC = point.cc pointbag.cc pntbag_i.cc

LIBOBJ = $(CXXSRC:%.cc=%.o)

INC =  point.h pointbag.h pointset.h pntbag_i.h

DEPENDINCLUDE = $(INC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libtopology

DISTFILES = $(CXXSRC) $(INC) Makefile LIBS.h

default:: $(DEPENDINCLUDE)
default:: $(TOPDIR)/lib/libtopology.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libtopology.so
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif
