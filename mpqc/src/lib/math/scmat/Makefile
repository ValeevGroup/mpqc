TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
DEFINES += -DSRCDIR=\"$(SRCDIR)\"

default:: $(DEPENDINCLUDE)
default:: $(TOPDIR)/lib/libscmat.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libscmat.so
endif

CXXSRC = \
         repl.cc       replrect.cc   repldiag.cc   replsymm.cc   replvect.cc \
         local.cc      localrect.cc  localdiag.cc  localsymm.cc  localvect.cc \
         abstract.cc   matrix.cc     matrix_i.cc \
         block.cc      blkiter.cc    elemop.cc     result.cc \
         matrix3.cc    vector3.cc    vector3_i.cc

CSRC = cmatrix.c svd.c

LIBOBJ = $(CXXSRC:%.cc=%.o) $(CSRC:%.c=%.o) 

INC =  cmatrix.h matrix_i.h cmatrix.h matrix.h abstract.h blkiter.h local.h \
       result.h vector3.h vector3_i.h block.h matrix3.h elemop.h repl.h

GENINC =

DEPENDINCLUDE = $(INC) $(GENINC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libscmat

TESTSRC = localtest.cc repltest.cc matrixtest.cc
TESTOBJ = matrixtest.o
TESTFILES = # localtest.in

DISTFILES = $(CXXSRC) $(INC) scmat.h $(WEBSRC) Makefile LIBS.h \
            $(TESTSRC) $(TESTFILES)

matrixtest:
	@echo The target must be a specialization such as localtest.

localtest: localtest.o $(TESTOBJ) \
	 $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(CXX) $(CXXFLAGS) -o localtest $^ $(SYSLIBS) 

repltest: repltest.o $(TESTOBJ) \
	$(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(CXX) $(CXXFLAGS) -o repltest $^ $(SYSLIBS) 

repltest.o: repltest.cc
	$(CXX) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $< -o $@

TAGS: $(CXXSRC) $(CSRC)
	etags $^

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTSRC:%.cc=%.d)
endif

