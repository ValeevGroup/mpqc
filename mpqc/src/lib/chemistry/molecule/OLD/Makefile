TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

CXXFLAGS = $(CXXFLAGS_DEBUG)

CXXSRC = atomcent.cc chemelem.cc molecule.cc \
         simple.cc symm.cc stre.cc bend.cc out.cc tors.cc lin.cc \
         simpleQCList.cc symmQCList.cc energy.cc coor.cc \
         shape.cc \

LIBOBJ = $(CXXSRC:%.cc=%.o)

INC = atomcent.h atomcentdefs.h chemelem.h molecule.h \
      simple.h symm.h localdef.h simpleQCList.h symmQCList.h \
      coor.h energy.h shape.h

DEPENDINCLUDE = $(INC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libmolecule

DISTFILES = $(CXXSRC) $(INC) Makefile

default:: $(DEPENDINCLUDE)
default:: libmolecule.a

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif
