TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

BIN_OR_LIB = LIB
TARGET_TO_MAKE = liboptions

LIBOBJ = GetLongOpt.o

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

sample: README.o liboptions.a
	$(CXX) -o $@ $^

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif


