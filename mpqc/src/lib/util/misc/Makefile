TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libmisc
BIN_OR_LIB = LIB

WEBFILES = pathfopen.web rdmatrix.web
.PRECIOUS: $(WEBFILES:%.web=%.c)

CSRC = ieee.c assert.c checkalloc.c mtype.c
GENSRC = $(WEBFILES:%.web=%.c)

CXXSRC = compute.cc error.cc formio.cc bug.cc regtime.cc

TESTSRC = comptest.cc formiot.cc bugtest.cc tregtime.cc

INC = assert.h libmisc.h memory.h timer.h compute.h ieee.h bug.h
GENINC = assert.gbl

DISTFILES = $(WEBFILES) $(CSRC) $(CXXSRC) $(TESTSRC) $(INC) Makefile LIBS.h

LIBOBJ= $(CSRC:%.c=%.$(OBJSUF)) \
        $(CXXSRC:%.cc=%.$(OBJSUF)) \
        $(GENSRC:%.c=%.$(OBJSUF))

TESTOBJ= $(TESTSRC:%.cc=%.$(OBJSUF))

DEPENDINCLUDE = $(INC) $(GENINC)

#####################################################################

default:: $(DEPENDINCLUDE)
default:: $(TOPDIR)/lib/libmisc.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libmisc.so
endif

comptest: comptest.cc $(TOPDIR)/lib/libmisc.a $(TOPDIR)/lib/libstate.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o comptest

tregtime: tregtime.cc $(TOPDIR)/lib/libmisc.a $(TOPDIR)/lib/libref.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o tregtime

bugtest: bugtest.cc $(TOPDIR)/lib/libmisc.a \
	$(TOPDIR)/lib/libkeyval.a $(TOPDIR)/lib/libstate.a \
	$(TOPDIR)/lib/libclass.a \
	$(TOPDIR)/lib/libcontainer.a $(TOPDIR)/lib/libref.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" $^ $(SYSLIBS) -o bugtest

formiot: formiot.cc libmisc.a
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ -o formiot

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d) $(TESTOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTOBJ:.o=.d)
endif

#########################################################################

ifneq ($(SRCDIR),.)
%.c: %.web
	cp $< .
	$(CTANGLE) $(notdir $<)
	rm -f $(notdir $<)

distclean::
	rm -f $(WEBFILES)
else
%.c: %.web
	$(CTANGLE) $<
endif

