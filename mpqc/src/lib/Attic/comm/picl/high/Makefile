TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = piclhigh
BIN_OR_LIB = LIB

COMMONSRC = piclgopf2.c gray.c f77stubs.c
DISTFILES = Makefile $(COMMONSRC) piclhcube.c piclparagon.c

LIBSRC = $(COMMONSRC)
PICLHIGH = none
ifeq ($(ARCH),PARAGON)
  PICLHIGH = piclparagon.c
endif
ifeq ($(PICLHIGH),none)
  PICLHIGH = piclhcube.c
endif

LIBSRC := $(LIBSRC) $(PICLHIGH)

LIBOBJ = $(LIBSRC:%.c=%.$(OBJSUF))

default:: ./piclhigh.a

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

