TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = piclext
BIN_OR_LIB = LIB

DISTFILES = Makefile gop.c prober.c loop.c ipscpicl.c gdcomb.c piclext.h

LIBSRC = gop.c prober.c loop.c ipscpicl.c
ifeq ($(ARCH),I860)
  LIBSRC := $(LIBSRC) gdcomb.c
endif

DEPENDCSRC = $(LIBSRC)
SABER_FILES = $(LIBSRC)

LIBOBJ = $(LIBSRC:%.c=%.$(OBJSUF))

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

default:: ./piclext.a

install_inc::
	$(INCINSTALL) piclext.h $(INCDIR)

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

