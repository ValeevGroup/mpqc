TOPDIR = ../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libmisc
BIN_OR_LIB = LIB

#CFLAGS := $(CFLAGS_DEBUG)
#FFLAGS := $(FFLAGS_DEBUG)

WEBFILES = debug.web pathfopen.web rdmatrix.web timer.web

CSRC = ieee.c prober.c assert.c machtype.c checkalloc.c mtype.c
GENSRC = $(WEBFILES:%.web=%.c)

CXXSRC = compute.cc compute_i.cc error.cc scostream.cc identity.cc

TESTSRC = comptest.cc scotest.cc

INC = assert.h libmisc.h memory.h timer.h compute.h scostream.h ieee.h
GENINC = machtype.gbl timer.gbl assert.gbl

DISTFILES = $(WEBFILES) $(CSRC) $(CXXSRC) $(TESTSRC) $(INC) Makefile LIBS.h

LIBOBJ= $(CSRC:%.c=%.$(OBJSUF)) \
        $(CXXSRC:%.cc=%.$(OBJSUF)) \
        $(GENSRC:%.c=%.$(OBJSUF))

DEPENDINCLUDE = $(INC) $(GENINC)

#####################################################################

default:: $(DEPENDINCLUDE)
default:: libmisc.a

comptest: comptest.cc libmisc.a ../state/libstate.a
	$(CXX) $(CXXFLAGS) $^ -o comptest

scotest: scotest.cc libmisc.a
	$(CXX) $(CXXFLAGS) $^ -o scotest

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

#########################################################################

ifneq ($(SRCDIR),.)
%.c: %.web
	cp $< .
	$(CTANGLE) $(notdir $<)
	rm -f $(notdir $<)
else
%.c: %.web
	$(CTANGLE) $<
endif

