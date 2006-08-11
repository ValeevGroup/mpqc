TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libSCoptions

LIBOBJ = GetLongOpt.$(OBJSUF)

TESTPROGS = sample

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

sample: README.$(OBJSUF) libSCoptions.$(LIBSUF)
	$(LTLINK) $(CXX) -o $@ $^ $(LTLINKBINOPTS)

$(LIBOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.$(OBJSUF)=.d)
endif


