
TOPDIR = ../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

DEFINES += -DSRCDIR=\"$(SRCDIR)\"

TARGET_TO_MAKE = libcints
BIN_OR_LIB = LIB

SRCS = nucrep.cc integraljf.cc ovlp.cc kei.cc pei.cc \
       fjtable.cc

LIBOBJ= $(SRCS:%.cc=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)
default:: $(TOPDIR)/lib/libcints.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libcints.so
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

ctest: ctest.o libcints.a $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LD) $(LDFLAGS) -o ctest $^ $(SYSLIBS)

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) ctest.d
endif

