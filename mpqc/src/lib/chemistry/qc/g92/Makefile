TOPDIR = ../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libg92
BIN_OR_LIB = LIB

SRCS = g92.cc g92scf.cc g92uhf.cc

CLASSES = Gaussian92 Gaussian92SCF

TESTSRCS = g92test.cc

LIBOBJ= $(SRCS:%.cc=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)
default:: dlinfo
default:: $(TOPDIR)/lib/libg92.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libg92.so
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

LIBS = $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

g92test:: g92test.o $(LIBS)
	$(CXX) $(LDFLAGS) -o g92test $^ $(SYSLIBS)

g92test.o: g92test.cc
	$(CXX) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTSRCS:%.cc=%.d)
endif

