TOPDIR=../../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libscf
BIN_OR_LIB = LIB

SRCS = fock.cc

CLASSES = Fock

LIBOBJ= $(SRCS:%.cc=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)
default:: dlinfo
default:: $(TOPDIR)/lib/libscf.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libscf.so
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

libs: $(SRCDIR)/LIBS.h
	$(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h > libs

scftest:: libs
	@echo "making scftest"

scftest:: scftest.o $(shell cat libs)
	$(CXX) $(LDFLAGS) -o scftest $^ $(SYSLIBS)

scftest.o: scftest.cc
	$(CXX) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) scftest.d
endif

