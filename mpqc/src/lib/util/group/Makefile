TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

REQUIREDCXXSRC = mstate.cc \
         message.cc messproc.cc messint.cc \
         reduce.cc topology.cc hcube.cc picl.cc
# pool.cc
# memory.cc

OPTIONALCXXSRC = messshm.cc messpgon.cc

CXXSRC := $(REQUIREDCXXSRC)
ifneq ($(ARCH),PARAGON)
  CXXSRC += messshm.cc
endif
ifeq ($(ARCH),PARAGON)
  CXXSRC += messpgon.cc
endif

INC = pool.h message.h mstate.h topology.h picl.h messshm.h
# memory.h

LIBOBJ = $(CXXSRC:%.cc=%.o)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libgroup

TESTCXXSRC = pooltest.cc messtest.cc
DISTFILES = $(REQUIREDCXXSRC) $(OPTIONALCXXSRC) $(INC) \
            Makefile LIBS.h $(TESTCXXSRC)

DEPENDINCLUDE =

default:: $(DEPENDINCLUDE)
default:: libgroup.a

pooltest: pooltest.o libgroup.a
	$(CXX) $(CXXFLAGS) -o pooltest $^

messtest.o: messtest.cc
	$(CXX) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

messtest: messtest.o \
	libgroup.a \
	libkeyval.a \
	libstate.a \
	libclass.a \
	libcontainer.a \
	libmisc.a \
	-lfl
	$(LD) $(LDFLAGS) -o messtest $^ $(SYSLIBS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

LD = $(CXX)

$(TESTCXXSRC:%.cc=%.d) $(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(TESTCXXSRC:%.cc=%.d) $(LIBOBJ:.o=.d)
endif

