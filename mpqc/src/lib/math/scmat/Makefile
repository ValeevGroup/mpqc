TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
DEFINES += -DSRCDIR=\"$(SRCDIR)\"

default:: $(DEPENDINCLUDE)

CXXSRC = dim.cc \
         abstract.cc   matrix.cc     matrix_i.cc \
         block.cc      blkiter.cc    elemop.cc     result.cc \
         local.cc      localrect.cc  localdiag.cc  localsymm.cc  localvect.cc \
         disthql.cc \
         dist.cc       distrect.cc   distdiag.cc   distsymm.cc   distvect.cc \
         repl.cc       replrect.cc   repldiag.cc   replsymm.cc   replvect.cc \
         blocked.cc \
         blockedvect.cc blockedrect.cc blockeddiag.cc blockedsymm.cc \
         matrix3.cc    vector3.cc    vector3_i.cc

TESTSRC = disttest.cc localtest.cc repltest.cc blockedtest.cc

CSRC = cmatrix.c svd.c

LIBOBJ = $(CXXSRC:%.cc=%.o) $(CSRC:%.c=%.o) 

INC =  cmatrix.h matrix_i.h cmatrix.h matrix.h abstract.h blkiter.h local.h \
       result.h vector3.h vector3_i.h block.h matrix3.h elemop.h repl.h \
       blocked.h disthql.h dim.h

GENINC =

DEPENDINCLUDE = $(INC) $(GENINC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libscmat

TESTOBJ = matrixtest.o

DISTFILES = $(CXXSRC) $(INC) scmat.h $(WEBSRC) Makefile LIBS.h \
            $(TESTSRC) $(TESTFILES)

LIBS = $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

matrixtest:
	@echo The target must be a specialization such as localtest.

localtest: localtest.o $(TESTOBJ) $(LIBS)
	$(CXX) $(CXXFLAGS) -o localtest $^ $(SYSLIBS) 

repltest: repltest.o $(TESTOBJ) $(LIBS)
	$(CXX) $(CXXFLAGS) -o repltest $^ $(SYSLIBS) 

disttest: disttest.o $(TESTOBJ) $(LIBS)
	$(CXX) $(CXXFLAGS) -o disttest $^ $(SYSLIBS) 

blockedtest: blockedtest.o $(TESTOBJ) $(LIBS)
	$(CXX) $(CXXFLAGS) -o blockedtest $^ $(SYSLIBS) 

repltest.o: repltest.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $< -o $@

disttest.o: disttest.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $< -o $@

TAGS: $(CXXSRC) $(CSRC)
	etags $^

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTSRC:%.cc=%.d)
endif

