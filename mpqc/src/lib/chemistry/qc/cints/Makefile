
TOPDIR=../../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

DEFINES += -DSRCDIR=\"$(SRCDIR)\"

TARGET_TO_MAKE = libcints
BIN_OR_LIB = LIB

SRCS = nucrep.cc integraljf.cc ovlp.cc kei.cc pei.cc \
       fjtable.cc int2jf.cc hrr.cc build.cc build_a0b0.cc top_build_a0b0.cc

LIBOBJ= $(SRCS:%.cc=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

ctest: ctest.o libcints.a $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LD) $(LDFLAGS) -o ctest $^ $(SYSLIBS)

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) ctest.d
endif

