TOPDIR = ../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

#CXXFLAGS = $(CXXFLAGS_DEBUG)

TARGET_TO_MAKE = libmpqc
BIN_OR_LIB = LIB

SRCS = mpqc.cc

TESTSRCS = mpqctest.cc

LIBOBJ= $(SRCS:%.cc=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)
default:: libmpqc.a

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

mpqctest: mpqctest.o libmpqc.a \
	$(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(CXX) $(LDFLAGS) -o mpqctest $^ $(SYSLIBS)

mpqctest.o: mpqctest.cc
	$(CXX) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTSRCS:%.cc=%.d)
endif

