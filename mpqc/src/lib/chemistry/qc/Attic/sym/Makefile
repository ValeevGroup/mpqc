TOPDIR = ../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

SGENFLAGS = -ifndef _chemistry_qc_sym_sym_h

TARGET_TO_MAKE = libsym
BIN_OR_LIB = LIB

CSRC = alloc.c create_r.c skeleton.c syminit.c aotoso.c mkrpd.c mksym.c \
       symerr.c sympack.c symutil.c

SGENFILE = symm.sgen
SGENSRC = symmallc.c symmasgn.c symmfree.c symmbc0.c symmrbc0.c symmsbc0.c \
            symmsnd0.c symmrcv0.c symmbrd.c symmbwr.c symmfrd.c symmfwr.c \
            symminit.c symmip.c symmiseq.c symmprnt.c symmzero.c

INC = sym_lib.h symm_mac.h 
SGENINC = $(SGENSRC:%.c=%.h)
GENINC = $(CSRC:%.c=%.gbl) $(CSRC:%.c=%.lcl)

LIBSRC = $(CSRC) $(SGENSRC)
LIBOBJ = $(LIBSRC:%.c=%.$(OBJSUF))

DISTFILES = $(CSRC) $(INC) symm.sgen Makefile LIBS.h

DEPENDINCLUDE = $(INC) $(SGENINC) $(GENINC)

###################################################################3

default:: $(DEPENDINCLUDE)
default:: libsym.a

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

