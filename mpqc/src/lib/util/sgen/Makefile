TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libsgen
BIN_OR_LIB = LIB

CSRC = chkma.c send0.c recv0.c sndrcv0.c iseq.c bcast0.c \
       rbcast0.c sbcast0.c print.c

INC = bcast0.h sndrcv0.h sgen.h

LIBOBJ = $(CSRC:%.c=%.$(OBJSUF))

DISTFILES = $(CSRC) $(INC) Makefile LIBS.h

##############################################

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

