TOPDIR = ../../../..
ifndef SRCDIR
  SRCDIR=.
endif

BASISDIR = $(shell $(SRCDIR)/$(TOPDIR)/bin/rmdotdot $(SRCDIR)/$(TOPDIR)/lib)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

DEFINES += -DBASISDIR=\"$(BASISDIR)\" -DSRCDIR=\"$(SRCDIR)\"

TARGET_TO_MAKE = libbasis
BIN_OR_LIB = LIB

CXXSRC = gaussbas.cc gaussshell.cc gaussshval.cc gaussbaval.cc \
       gaussbas_i.cc gaussshe_i.cc

INC = gaussbas.h gaussshe_i.h basis.h gaussbas_i.h gaussshell.h

TESTSRC = btest.cc
TESTOBJ = $(TESTSRC:%.cc=%.o)
TESTFILES = btest.in

DISTFILES = $(CXXSRC) $(INC) Makefile $(TESTSRC) $(TESTFILES)

LIBOBJ= $(CXXSRC:%.cc=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)
default:: libbasis.a

btest: $(TESTOBJ) libbasis.a \
	 $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(CXX) $(CXXFLAGS) -o btest $^ $(SYSLIBS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTOBJ:%.o=%.d)
endif

