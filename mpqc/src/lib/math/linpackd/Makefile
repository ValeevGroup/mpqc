TOPDIR = ../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = linpackd
BIN_OR_LIB = LIB

CSRC = linpackd.c

LIBSRC = $(CSRC)

DISTFILES = $(CSRC) Makefile LIBS.h

LIBOBJ := $(LIBSRC:%.c=%.$(OBJSUF))

default:: $(TOPDIR)/lib/linpackd.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/linpackd.so
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

DFILES := $(LIBSRC:.c=.d)
ifneq ($(DODEPEND),no)
ifneq ($(DFILES),)
include $(DFILES)
endif
endif

