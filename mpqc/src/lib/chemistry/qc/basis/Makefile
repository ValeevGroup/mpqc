TOPDIR = ../../../..
ifndef SRCDIR
  SRCDIR=.
endif

BASISDIR = $(shell $(SRCDIR)/$(TOPDIR)/bin/rmdotdot $(SRCDIR)/$(TOPDIR)/lib)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

DEFINES += -DBASISDIR=\"$(BASISDIR)/basis\" -DSRCDIR=\"$(SRCDIR)\"

TARGET_TO_MAKE = libbasis
BIN_OR_LIB = LIB

CXXSRC = gaussbas.cc gaussshell.cc gaussshval.cc gaussbaval.cc \
       gaussbas_i.cc gaussshe_i.cc files.cc

INC = gaussbas.h gaussshe_i.h basis.h gaussbas_i.h gaussshell.h files.h

TESTSRC = btest.cc
TESTOBJ = $(TESTSRC:%.cc=%.o)
TESTFILES = btest.kv

DISTFILES = $(CXXSRC) $(INC) Makefile $(TESTSRC) $(TESTFILES)

LIBOBJ= $(CXXSRC:%.cc=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)
default:: $(TOPDIR)/lib/libbasis.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libbasis.so
endif

btest: $(TESTOBJ) $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LD) $(LDFLAGS) -o btest $^ $(SYSLIBS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTOBJ:%.o=%.d)
endif

