TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

LD=$(CXX)

CXXSRC = nlp.cc opt.cc update.cc qnewton.cc gdiis.cc efc.cc dfp.cc powell.cc \
         scextrap.cc scextrapmat.cc diis.cc
CSRC = 

CLASSES = Optimize QNewtonOpt GDIISOpt LineOpt EFCOpt \
          NLP0 NLP1 NLP2 \
          HessianUpdate DFPUpdate BFGSUpdate PowellUpdate

LIBOBJ = $(CXXSRC:%.cc=%.o) $(CSRC:%.c=%.o)

INC = nlp.h opt.h gdiis.h efc.h update.h qnewton.h \
      scextrap.h scextrapmat.h diis.h

DEPENDINCLUDE =

BIN_OR_LIB = LIB
TARGET_TO_MAKE = liboptimize

DISTFILES = opttest.cc $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h

default:: $(DEPENDINCLUDE)
default:: dlinfo
default:: $(TOPDIR)/lib/liboptimize.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/liboptimize.so
endif

opttest: opttest.o optimize.a util/options/liboptions.a \
	$(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LD) $(LDFLAGS) -o opttest $^ $(SYSLIBS)

opttest.o: opttest.cc
	$(CXX) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

scextest: scextest.o liboptimize.a \
	$(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LD) $(LDFLAGS) -o scextest $^ $(SYSLIBS)

scextest.o: scextest.cc
	$(CXX) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

kvopt: kvopt.o $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LD) $(LDFLAGS) -o kvopt $^ $(SYSLIBS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

kvopt.d opttest.d scextest.d $(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) kvopt.d opttest.d scextest.d
endif

