TOPDIR=../..
ifndef SRCDIR
  SRCDIR=.
endif

BIN_OR_LIB = BIN

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

LD = $(CXX)

TARGET_TO_MAKE = molrender

CXXSRC = main.cc
CSRC =
INC =

DEPENDINCLUDE = $(INC)

BINOBJ = $(CXXSRC:%.cc=%.o) $(CSRC:%.c=%.o)

DISTFILES = $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h molrender.in

default:: molrender

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

main.o: main.cc
	$(CXX) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

$(BINOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.o=.d)
endif

