TOPDIR=../../../../..
ifndef SRCDIR
  SRCDIR=.
endif

BASISDIR = $(shell $(SRCDIR)/$(TOPDIR)/bin/rmdotdot $(SRCDIR)/$(TOPDIR)/lib)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

DEFINES += -DBASISDIR=\"$(BASISDIR)/basis\" -DSRCDIR=\"$(SRCDIR)\"

TARGET_TO_MAKE = libbasis
BIN_OR_LIB = LIB

CXXSRC = gaussbas.cc gaussbas_i.cc gaussshell.cc gaussshe_i.cc files.cc \
         gaussbaval.cc gaussshval.cc cartiter.cc cartiter_i.cc \
         transform.cc shellrot.cc integral.cc obint.cc petite.cc aotoso.cc \
         symmint.cc ortho.cc tbint.cc dercent.cc

INC = basis.h gaussbas.h guassshell.h gaussbas_i.h gaussshe_i.h files.h \
      cartiter.h cartiter_i.h transform.h shellrot.h \
      integral.h obint.h petite.h symmint.h tbint.h dercent.h

TESTSRC = btest.cc
TESTOBJ = $(TESTSRC:%.cc=%.o)
TESTFILES = btest.kv

DISTFILES = $(CXXSRC) $(INC) Makefile $(TESTSRC) $(TESTFILES)

LIBOBJ= $(CXXSRC:%.cc=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)
default:: $(TOPDIR)/lib/libbasis.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libbasis.so
endif

btest: $(TESTOBJ) libintv3.a liboint3.a libintv2.a liboint2.a \
                  $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LD) $(LDFLAGS) -o btest $^ $(SYSLIBS)

clean::
	rm -f btest btest.out

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTOBJ:%.o=%.d)
endif

