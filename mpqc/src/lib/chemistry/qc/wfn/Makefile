TOPDIR=../../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libwfn
BIN_OR_LIB = LIB

SRCS = wfn.cc densval.cc natorbval.cc obwfn.cc density.cc hcorewfn.cc \
       accum.cc effh.cc hcore.cc

CLASSES = ElectronDensity OneBodyWavefunction Wavefunction

LIBOBJ= $(SRCS:%.cc=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)
default:: dlinfo
default:: $(TOPDIR)/lib/libwfn.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libwfn.so
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

libs: $(SRCDIR)/LIBS.h
	$(LISTLIBS) $(MPQCLIBSDEFINE) $(INCLUDE) $(SRCDIR)/LIBS.h > libs

wfntest:: libs
	@echo "making wfntest"

wfntest:: wfntest.o $(shell cat libs)
	$(CXX) $(LDFLAGS) -o wfntest $^ $(SYSLIBS)

wfntest.o: wfntest.cc
	$(CXX) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

