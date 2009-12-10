TOPDIR=../../../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
include $(TOPDIR)/lib/Makedirlist

TARGET_TO_MAKE = libSCpsi
BIN_OR_LIB = LIB

CXXSRC = psiwfn.cc psiexenv.cc psiinput.cc psifile11.cc psicc.cc psici.cc psiqtorder.cc \
r12ref.cc
CSRC = 
INC = psiwfn.h psiinput.h psiexenv.h psifile11.h psicc.h psici.h psiqtorder.h \
r12ref.h

# R12 functionality depends on MP2-R12 code
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_MBPTR12),yes)
CXXSRC := $(CXXSRC) psicc_pt2r12.cc rdm.cc
#psici_pt2r12.cc psiwfn_pt2r12.cc 
INC := $(INC) psicc_pt2r12.h
#psici_pt2r12.h
endif

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

