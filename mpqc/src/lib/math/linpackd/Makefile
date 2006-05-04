TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libSClinpackd
BIN_OR_LIB = LIB

CSRC = linpackd.c

LIBSRC = $(CSRC)

DISTFILES = $(CSRC) Makefile LIBS.h

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

LIBOBJ := $(LIBSRC:%.c=%.$(OBJSUF))

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

DFILES := $(LIBSRC:.c=.d)
ifneq ($(DODEPEND),no)
ifneq ($(DFILES),)
include $(DFILES)
endif
endif

