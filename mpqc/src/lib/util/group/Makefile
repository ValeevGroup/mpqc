TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

CXXSRC = pool.cc mstate.cc \
         message.cc messproc.cc messint.cc messshm.cc
# memory.cc

INC = pool.h message.h mstate.h
# memory.h

LIBOBJ = $(CXXSRC:%.cc=%.o)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libgroup

TESTCXXSRC = pooltest.cc messtest.cc
DISTFILES = $(CXXSRC) $(INC) Makefile LIBS.h $(TESTCXXSRC)

DEPENDINCLUDE =

default:: $(DEPENDINCLUDE)
default:: libgroup.a

pooltest: pooltest.o libgroup.a
	$(CXX) $(CXXFLAGS) -o pooltest $^

messtest: messtest.o \
	libgroup.a \
	../../util/keyval/libkeyval.a \
	../../util/state/libstate.a \
	../../util/class/libclass.a \
	../../util/container/libcontainer.a \
	../../util/misc/libmisc.a \
	-lfl
	$(CXX) $(CXXFLAGS) -o messtest $^

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(TESTCXXSRC:%.cc=%.d) $(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(TESTCXXSRC:%.cc=%.d) $(LIBOBJ:.o=.d)
endif

