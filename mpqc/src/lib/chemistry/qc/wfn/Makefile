TOPDIR=../../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libwfn
BIN_OR_LIB = LIB

SRCS = wfn.cc densval.cc natorbval.cc obwfn.cc density.cc hcorewfn.cc \
       accum.cc hcore.cc

CLASSES = ElectronDensity OneBodyWavefunction Wavefunction

LIBOBJ= $(SRCS:%.cc=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

LIBS = $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

wfntest:: wfntest.o $(LIBS)
	$(CXX) $(LDFLAGS) -o wfntest $^ $(SYSLIBS)

wfntest.o: wfntest.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

