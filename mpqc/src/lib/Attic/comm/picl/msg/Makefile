TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

BIN_OR_LIB = LIB
TARGET_TO_MAKE = piclmsg

LIBSRC = piclmsg.c
DISTFILES = Makefile $(LIBSRC)

LIBOBJ = $(LIBSRC:%.c=%.$(OBJSUF))

default:: ./piclmsg.a

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

