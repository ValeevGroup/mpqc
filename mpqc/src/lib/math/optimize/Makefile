TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

LD=$(CXX)

CXXSRC = function.cc opt.cc update.cc qnewton.cc \
         gdiis.cc efc.cc dfp.cc powell.cc \
         scextrap.cc scextrapmat.cc diis.cc transform.cc conv.cc
CSRC = 

CLASSES = Optimize QNewtonOpt GDIISOpt LineOpt EFCOpt \
          Function \
          HessianUpdate DFPUpdate BFGSUpdate PowellUpdate

LIBOBJ = $(CXXSRC:%.cc=%.o) $(CSRC:%.c=%.o)

INC = function.h opt.h gdiis.h efc.h update.h qnewton.h \
      scextrap.h scextrapmat.h diis.h transform.h conv.h

DEPENDINCLUDE =

BIN_OR_LIB = LIB
TARGET_TO_MAKE = liboptimize

DISTFILES = opttest.cc $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h

default:: $(DEPENDINCLUDE)

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

