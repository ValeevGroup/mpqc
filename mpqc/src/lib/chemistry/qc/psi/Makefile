TOPDIR=../../../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libSCpsi
BIN_OR_LIB = LIB

SRCS = psiwfn.cc psiexenv.cc psiinput.cc psifile11.cc
CSRCS = 

psitest.testrun: DO_TESTRUN=no

TESTSRCS = psitest.cc

TESTPROGS = psitest

LIBOBJ = $(SRCS:%.cc=%.$(OBJSUF))
LIBOBJ += $(CSRCS:%.c=%.$(OBJSUF))

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

