TOPDIR=../../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libpsi
BIN_OR_LIB = LIB

SRCS = psi.cc file11.cc psiinput.cc psici.cc psiccsd.cc psiccsd_t.cc psiccsdt.cc
CSRCS = read11.c

TESTSRCS = psitest.cc

LIBOBJ = $(SRCS:%.cc=%.$(OBJSUF))
LIBOBJ += $(CSRCS:%.c=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)
interface:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

psitest: psitest.o  \
	$(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(CXX) $(LDFLAGS) -o psitest $^ $(SYSLIBS)

psitest.o: psitest.cc
	$(CXX) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

print: 
	echo $(SRCDIR)

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTSRCS:%.cc=%.d)
endif

