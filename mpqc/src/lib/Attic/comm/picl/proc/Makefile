TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = piclproc
BIN_OR_LIB = LIB

LIBSRC = piclproc.c
DISTFILES = Makefile $(LIBSRC)

LIBOBJ = $(LIBSRC:%.c=%.$(OBJSUF))

default:: ./piclproc.a

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

