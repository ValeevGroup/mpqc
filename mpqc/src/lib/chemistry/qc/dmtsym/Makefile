TOPDIR=../../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

SGENFLAGS = -ifndef _chemistry_qc_dmtsym_sym_h

TARGET_TO_MAKE = libdmtsym
BIN_OR_LIB = LIB

CSRC = create_r.c skeleton.c syminit.c mkrpd.c mksym.c

CXXSRC = sym_cpp.cc

SGENFILE = symm.sgen
SGENSRC = symmallc.c symmasgn.c symmfree.c symmbc0.c symmrbc0.c symmsbc0.c \
          symmsnd0.c symmrcv0.c symminit.c symmiseq.c symmprnt.c symmzero.c

INC = sym_dmt.h symm_mac.h
SGENINC = $(SGENSRC:%.c=%.h)
GENINC = $(CSRC:%.c=%.gbl) $(CSRC:%.c=%.lcl)

LIBSRC = $(CSRC) $(SGENSRC) $(CXXSRC)
LIBOBJ = $(LIBSRC:%.c=%.$(OBJSUF))
LIBOBJ := $(LIBOBJ:%.cc=%.$(OBJSUF))

DISTFILES = $(CSRC) $(INC) symm.sgen Makefile LIBS.h

DEPENDINCLUDE = $(SGENINC) $(GENINC) $(INC)

#####################################################################

default:: $(DEPENDINCLUDE)
default:: $(TOPDIR)/lib/libdmtsym.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libdmtsym.so
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

