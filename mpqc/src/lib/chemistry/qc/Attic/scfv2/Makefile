TOPDIR = ../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

SGENFLAGS = -ifndef _chemistry_qc_scfv2_scf_h

TARGET_TO_MAKE = libscfv2
BIN_OR_LIB = LIB

ifeq ($(ARCH),SUN4)
override COPTIONS = -O2
endif

DEFINES += -DBASISLIB="\"$(LIBDIR)/basis.ipv2\""

CSRC = scf_core.c scf_diis.c scf_en.c scf_file.c scf_fock.c scf_gmat.c \
       scf_init.c scf_inp.c scf_iter.c scf_mkden.c scf_mkg.c scf_mkgd.c \
       scf_mkpk.c scf_oeis.c scf_orth.c scf_vect.c

SGENFILE = scf.sgen
SGENSRC = scfallc.c scfasgn.c scfbc0.c scfbrd.c scfbwr.c scffree.c scfinit.c \
          scfiseq.c scfprnt.c scfrbc0.c scfsbc0.c scfzero.c scfsnd0.c scfrcv0.c

INC = scf_lib.h
SGENINC = $(SGENSRC:%.c=%.h)
GENINC = $(CSRC:%.c=%.gbl) $(CSRC:%.c=%.lcl)


LIBSRC = $(CSRC) $(SGENSRC)
LIBOBJ = $(LIBSRC:%.c=%.$(OBJSUF))

DEPENDINCLUDE = $(INC) $(SGENINC) $(GENINC)

DISTFILES = $(CSRC) $(INC) scf.sgen Makefile LIBS.h

############################################################

default:: $(DEPENDINCLUDE)
default:: libscfv2.a

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

