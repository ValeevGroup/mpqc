TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
DEFINES += -DSRCDIR=\"$(SRCDIR)\" -DSRCLIBDIR=\"$(SRCLIBDIR)\"

CXXSRC = atomcent.cc chemelem.cc molecule.cc \
         simple.cc stre.cc bend.cc out.cc tors.cc stors.cc linip.cc linop.cc \
         energy.cc coor.cc imcoor.cc symmcoor.cc redund.cc cartcoor.cc \
         molshape.cc chmelm_i.cc taylor.cc \
         molrender.cc molinfo.cc

CLASSES = Molecule \
          AtomicCenter ChemicalElement \
          IntCoor SumIntCoor SetIntCoor \
          MolecularCoor IntMolecularCoor \
          MolecularEnergy  TaylorMolecularEnergy \
          VDWShape ConnollyShape \
          SimpleCo StreSimpleCo BendSimpleCo TorsSimpleCo \
          OutSimpleCo LinIPSimpleCo LinOPSimpleCo \
          MolInfo AtomInfo \
          RenderedMolecule RenderedStickMolecule RenderedBallMolecule \
          RenderedMolecularSurface

LIBOBJ = $(CXXSRC:%.cc=%.o)

INC = atomcent.h atomcentdefs.h chemelem.h molecule.h \
      simple.h localdef.h \
      coor.h energy.h molshape.h chmelm_i.h taylor.h taylor_f.h \
      molrender.h molinfo.h

DEPENDINCLUDE = $(INC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libmolecule

TESTSRC = moltest.cc
TESTOBJ = $(TESTSRC:%.cc=%.o)
TESTFILES = moltest.in

DISTFILES = $(CXXSRC) $(INC) Makefile $(TESTSRC) $(TESTFILES)

default:: $(DEPENDINCLUDE)
default:: dlinfo
default:: $(TOPDIR)/lib/libmolecule.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libmolecule.so
endif
default:: libs

interface:: $(DEPENDINCLUDE)

LIBS = $(shell cat libs)

libs:
	 $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h > libs

LD = $(CXX)

moltest: $(TESTOBJ) $(LIBS)
	$(LD) $(LDFLAGS) -o moltest $^ $(SYSLIBS)

symmetrize: symmetrize.o $(LIBS)
	$(LD) $(LDFLAGS) -o symmetrize $^ $(SYSLIBS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTOBJ:%.o=%.d)
endif
