TOPDIR=../../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libscf
BIN_OR_LIB = LIB

SRCS = scf.cc scfvector.cc scfgradient.cc \
       clscf.cc

CLASSES = SCF CLSCF

LIBOBJ= $(SRCS:%.cc=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)
default:: dlinfo
default:: $(TOPDIR)/lib/libscf.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libscf.so
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

LIBS = $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

scftest:: scftest.o $(LIBS)
	$(CXX) $(LDFLAGS) -o scftest $^ $(SYSLIBS)

scftest.o: scftest.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) scftest.d
endif

