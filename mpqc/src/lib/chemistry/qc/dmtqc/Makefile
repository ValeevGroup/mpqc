TOPDIR = ../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libdmtqc
BIN_OR_LIB = LIB

CSRC = dens.c opdens.c ewdensity.c expewdens.c
INC = libdmtqc.h

LIBOBJ= $(CSRC:%.c=%.$(OBJSUF))

GENINC = $(CSRC:%.c=%.gbl)

DEPENDINCLUDE = $(INC) $(GENINC)

DISTFILES = $(CSRC) $(INC) Makefile LIBS.h

#############################################################

default:: $(DEPENDINCLUDE)
default:: $(TOPDIR)/lib/libdmtqc.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libdmtqc.so
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

