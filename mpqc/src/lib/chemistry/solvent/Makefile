TOPDIR=../../../..
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

interface:: $(DEPENDINCLUDE)

LIBS = $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

bemtest: $(TESTOBJ) $(LIBS)
	$(LD) $(LDFLAGS) -o bemtest $^ $(SYSLIBS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTOBJ:%.o=%.d)
endif
