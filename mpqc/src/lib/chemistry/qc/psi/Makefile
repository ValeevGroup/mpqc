TOPDIR=../../../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libSCpsi
BIN_OR_LIB = LIB

CXXSRC = psiwfn.cc psiexenv.cc psiinput.cc psifile11.cc psicc.cc psicc_pt2r12.cc
CSRC = 
INC = psiwfn.h psiinput.h psiexenv.h psifile11.h psicc.h psicc_pt2r12.h

psitest.testrun: DO_TESTRUN=no

TESTSRC = psitest.cc
TESTPROGS = psitest

LIBSRC = $(CXXSRC) $(CSRC)
LIBOBJ = $(LIBSRC:%.c=%.$(OBJSUF))
LIBOBJ := $(LIBOBJ:%.cc=%.$(OBJSUF))

DISTFILES = $(TESTSRC) $(INC) Makefile LIBS.h
DEPENDINCLUDE = $(INC) $(GENINC) $(SGENINC)

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

default:: $(DEPENDINCLUDE)
interface:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

psitest: psitest.$(OBJSUF)  \
	$(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LTLINK) $(CXX) $(LDFLAGS) -o psitest $^ $(SYSLIBS) $(LTLINKBINOPTS)

psitest.$(OBJSUF): psitest.cc
	$(LTCOMP) $(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

print: 
	echo $(SRCDIR)

$(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d) $(TESTSRCS:%.cc=%.d)
endif

