TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
DEFINES += -DSRCDIR=\"$(SRCDIR)\"

LD = $(CXX)
CXXSRC = bem.cc
CLASSES = BEMSolvent

LIBOBJ = $(CXXSRC:%.cc=%.o)

INC = bem.h

DEPENDINCLUDE = $(INC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libsolvent

TESTSRC = bemtest.cc
TESTOBJ = $(TESTSRC:%.cc=%.o)
TESTFILES = bemtest.in

DISTFILES = $(CXXSRC) $(INC) Makefile $(TESTSRC) $(TESTFILES)

default:: $(DEPENDINCLUDE)
default:: dlinfo
default:: $(TOPDIR)/lib/libsolvent.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libsolvent.so
endif

interface:: $(DEPENDINCLUDE)

bemtest: $(TESTOBJ) $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LD) $(LDFLAGS) -o bemtest $^ $(SYSLIBS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTOBJ:%.o=%.d)
endif
