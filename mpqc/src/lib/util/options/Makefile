TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

BIN_OR_LIB = LIB
TARGET_TO_MAKE = liboptions

LIBOBJ = GetLongOpt.o

default:: $(TOPDIR)/lib/liboptions.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/liboptions.so
endif

sample: README.o liboptions.a
	$(CXX) -o $@ $^

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif


