TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = picln6400
BIN_OR_LIB = LIB

LIBSRC = picln6400.c
DISTFILES = Makefile $(LIBSRC)

LIBOBJ = $(LIBSRC:%.c=%.$(OBJSUF))

default:: ./picln6400.a

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

