#!/net/munin/usr/local/bin/gmake
TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

#COPTIONS = $(COPTIONS_DEBUG)
#CXXOPTIONS = $(CXXOPTIONS_DEBUG)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = nlp

default:: nlp.a

#
#
# Use the AT&T compiler
#
#CXX = CC -DATandT
#NEWMAT = newmat7/newmat7.a 
#
# Use the gnu complier
#
#CXX = /net/munin/usr/local/bin/g++ -DGNUG
#NEWMAT = /usr/local/lib/newmat7.a /usr/local/lib/libg++.a
#
#

NEWTOBJS = optnewton.o newton.o qnewton.o \
           update.o upbfgs.o upfd.o upno.o
OBJS  = tstfcn.o nlp.o optimize.o optcg.o nlcg.o deriv.o search.o

LIBOBJ = $(OBJS) $(NEWTOBJS)

LIBS = $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

TESTSRC = opttest.cc testit.cc
TESTFILES = opttest.in

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include tstopt.d $(LIBOBJ:.o=.d) $(TESTSRC:%.cc=%.d)
endif

#
#  Main test routine
#

targetclean::
	rm -f topt tlnsrch tvec tform testit

topt:   tstopt.o nlp.a $(LIBS)
	$(CXX) $(DEFINES) -o topt $^

#
# Miscellaneous test routines
#

opttest: opttest.o nlp.a \
	 $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(CXX) -o opttest $^

tlnsrch: tstsrch.o $(LIBS)
	$(CXX) -o tlnsrch $^

tvec: tstvec.o $(LIBS)
	$(CXX) -o tvec $^

tform: tstform.o $(LIBS)
	$(CXX) -o tform $^

testit: testit.o ljeval.o
	$(CC) -o testit $^ -lm

lint: $(SRCS)
	lint $(SRCS)

