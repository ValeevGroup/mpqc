TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

CXXSRC = dvector3.cc dvector.cc dmatrix.cc lmath.cc newmat.cc

LIBOBJ = $(CXXSRC:%.cc=%.o)

INC =  nihmatrix.h lmath.h

DEPENDINCLUDE = $(INC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libnihmatrix

TESTSRC = matrixtest.cc
TESTOBJ = $(TESTSRC:%.cc=%.o)
TESTFILES = matrixtest.in

DISTFILES = $(CXXSRC) $(INC) Makefile LIBS.h $(TESTSRC) $(TESTFILES)

default:: $(DEPENDINCLUDE)
default:: libnihmatrix.a

matrixtest: $(TESTOBJ) libnihmatrix.a \
	 $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(CXX) $(CXXFLAGS) -o matrixtest $^

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTOBJ:%.o=%.d)
endif

