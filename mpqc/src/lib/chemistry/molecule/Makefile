TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
DEFINES += -DSRCDIR=\"$(SRCDIR)\"

CXXSRC = atomcent.cc chemelem.cc molecule.cc \
         simple.cc stre.cc bend.cc out.cc tors.cc lin.cc \
         energy.cc coor.cc \
         shape.cc chmelm_i.cc imcoor.cc taylor.cc

LIBOBJ = $(CXXSRC:%.cc=%.o)

INC = atomcent.h atomcentdefs.h chemelem.h molecule.h \
      simple.h localdef.h \
      coor.h energy.h shape.h chmelm_i.h taylor.h taylor_f.h

DEPENDINCLUDE = $(INC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libmolecule

TESTSRC = moltest.cc
TESTOBJ = $(TESTSRC:%.cc=%.o)
TESTFILES = moltest.in

DISTFILES = $(CXXSRC) $(INC) Makefile $(TESTSRC) $(TESTFILES)

default:: $(DEPENDINCLUDE)
default:: libmolecule.a

interface:: $(DEPENDINCLUDE)

moltest: $(TESTOBJ) libmolecule.a \
	 $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(CXX) $(CXXFLAGS) -o moltest $^ $(SYSLIBS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTOBJ:%.o=%.d)
endif
