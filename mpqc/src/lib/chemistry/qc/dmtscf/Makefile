TOPDIR = ../../../..
ifndef SRCDIR
  SRCDIR=.
endif

SGENFLAGS = -ifndef _chemistry_qc_dmtscf_scf_h

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libdmtscf
BIN_OR_LIB = LIB

CSRC = scf_bnd.c scf_core.c scf_diis.c scf_dip.c scf_en.c scf_esp.c \
       scf_fock.c scf_gmat.c scf_init.c scf_inp.c scf_iter.c scf_loopg.c \
       scf_lowd.c scf_mkden.c scf_mkgd.c scf_mull.c scf_oeis.c scf_orth.c \
       scf_proj.c scf_vect.c scf_ex.c scf_mkgdlb.c

SGENFILE = scf.sgen

SGENSRC = scfasgn.c scfbc0.c scfbrd.c scfbwr.c scfinit.c \
          scfiseq.c scfprnt.c scfrbc0.c scfsbc0.c scfsnd0.c scfrcv0.c \
          scfallc.c scffree.c scfzero.c
GENSRC = $(SGENSRC)

INC = scf_dmt.h masses.h vdw_rad.h

SGENINC = $(GENSRC:%.c=%.h)
GENINC = $(CSRC:%.c=%.gbl) $(CSRC:%.c=%.lcl)

LIBOBJ = $(CSRC:%.c=%.o) $(GENSRC:%.c=%.o)

DEPENDINCLUDE = $(INC) $(GENINC) $(SGENINC)

DISTFILES = $(CSRC) $(INC) scf.sgen Makefile LIBS.h

###################################################################

default:: $(DEPENDINCLUDE)
default:: libdmtscf.a

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

###################################################################

# The SGI cc compiler cannot handle scf_mkgd.c.
# The resulting code will sometimes produce a wrong answer.
ifeq ($(ARCH),SGI)
scf_mkgd.o: scf_mkgd.c
	gcc -O -g $(INCLUDE) -c $< -o scf_mkgd.o
scf_ex.o: scf_ex.c
	gcc -g $(INCLUDE) -c $< -o scf_ex.o
endif
