TOPDIR = ../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = linpackd
BIN_OR_LIB = LIB

FSRC = linpackd.f
CSRC = upper.c lower.c

LIBSRC = $(FSRC)

ifeq ($(MAP_C_F),upper)
  LIBSRC += upper.c
endif
ifeq ($(MAP_C_F),lower)
  LIBSRC += lower.c
endif

DISTFILES = $(FSRC) $(CSRC) Makefile LIBS.h

ifeq ($(FC),f2c)
GENSRC = $(FSRC:%.f=%.c) $(FSRC:%.f=%.c.c)
endif

LIBOBJ = $(LIBSRC:%.f=%.$(OBJSUF))
LIBOBJ := $(LIBOBJ:%.c=%.$(OBJSUF))

default:: linpackd.a

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

DFILES := $(LIBSRC:%.f=)
DFILES := $(DFILES:.c=.d)
ifneq ($(DODEPEND),no)
ifneq ($(DFILES),)
include $(DFILES)
endif
endif

