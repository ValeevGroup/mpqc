TOPDIR = ../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libbio
BIN_OR_LIB = LIB

ifeq ($(ARCH),RS6000)
DEFINES += -D_POSIX_SOURCE
endif

CSRC = bioinit.c bioopen.c errors.c file_info.c r_async.c ram.c s_async.c \
       sequential.c

LIBOBJ = $(CSRC:%.c=%.$(OBJSUF))

INC = libbio.h param.h pointers.h types.h
GENINC = $(CSRC:%.c=%.gbl)

DEPENDINCLUDE = $(INC) $(GENINC)

DISTFILES = $(CSRC) $(INC) Makefile LIBS.h

####################################################

default:: $(DEPENDINCLUDE)
default:: libbio.a

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

