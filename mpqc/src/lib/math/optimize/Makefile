TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

LD=$(CXX)

CXXSRC = nlp.cc opt.cc update.cc qnewton.cc gdiis.cc
CSRC = 

LIBOBJ = $(CXXSRC:%.cc=%.o) $(CSRC:%.c=%.o)

INC = nlp.h opt.h gdiis.h

DEPENDINCLUDE =

BIN_OR_LIB = LIB
TARGET_TO_MAKE = optimize

DISTFILES = opttest.cc $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h

default:: $(DEPENDINCLUDE)
default:: $(TOPDIR)/lib/optimize.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/optimize.so
endif

opttest: opttest.o optimize.a util/options/liboptions.a \
	$(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LD) $(LDFLAGS) -o opttest $^ $(SYSLIBS)

opttest.o: opttest.cc
	$(CXX) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

opttest.d $(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) opttest.d
endif

