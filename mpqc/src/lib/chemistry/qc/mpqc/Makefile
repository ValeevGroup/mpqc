TOPDIR=../../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libmpqc
BIN_OR_LIB = LIB

SRCS = mpqc.cc

CLASSES = MPSCF

TESTSRCS = mpqctest.cc

LIBOBJ= $(SRCS:%.cc=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)
default:: dlinfo
default:: $(TOPDIR)/lib/libmpqc.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libmpqc.so
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

LIBS = $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

mpqctest: mpqctest.o $(LIBS)
	$(CXX) $(LDFLAGS) -o mpqctest $^ $(SYSLIBS)

mpqctest.o: mpqctest.cc
	$(CXX) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTSRCS:%.cc=%.d)
endif

