TOPDIR = ../../../..
ifndef SRCDIR
  SRCDIR=.
endif

SGENFLAGS = -ifndef _chemistry_qc_dmtscf_scf_h

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libdmtscf
BIN_OR_LIB = LIB

CSRC = scf_bnd.c scf_core.c scf_diis.c scf_en.c \
       scf_gmat.c scf_iter.c scf_loopg.c \
       scf_mkden.c scf_mkgd.c scf_oeis.c scf_orth.c \
       scf_proj.c scf_vect.c scf_ex.c scf_mkgdlb.c scf_mkk.c scf_mkj.c \
       scf_loopj.c scf_loopk.c

CXXSRC = scf_kinp.cc

SGENFILE = scf.sgen

SGENSRC = scfasgn.c scfbc0.c scfinit.c \
          scfiseq.c scfprnt.c scfrbc0.c scfsbc0.c scfsnd0.c scfrcv0.c \
          scfallc.c scffree.c scfzero.c
GENSRC = $(SGENSRC)

INC = scf_dmt.h

SGENINC = $(GENSRC:%.c=%.h)
GENINC = $(CSRC:%.c=%.gbl) $(CSRC:%.c=%.lcl)

LIBOBJ = $(CSRC:%.c=%.o) $(CXXSRC:%.cc=%.o) $(GENSRC:%.c=%.o)

DEPENDINCLUDE = $(INC) $(GENINC) $(SGENINC)

DISTFILES = $(CSRC) $(INC) scf.sgen Makefile LIBS.h

###################################################################

default:: $(DEPENDINCLUDE)
default:: $(TOPDIR)/lib/libdmtscf.a
interface:: $(DEPENDINCLUDE)
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libdmtscf.so
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

ifeq ($(ARCH),PARAGON)
  MPQCLIBSDEFINE=-DUSE_PICLPARAGON
endif

dmtscftest: dmtscftest.o libdmtscf.a \
        $(shell $(LISTLIBS) $(MPQCLIBSDEFINE) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(CXX) $(LDFLAGS) -o $@ $^ $(SYSLIBS)

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) dmtscftest.d
endif

###################################################################

# The SGI cc compiler cannot handle scf_mkgd.c.
# The resulting code will sometimes produce a wrong answer.
#ifeq ($(ARCH),SGI)
#scf_mkgd.o: scf_mkgd.c
#	gcc -O -g $(INCLUDE) -c $< -o scf_mkgd.o
#scf_ex.o: scf_ex.c
#	gcc -g $(INCLUDE) -c $< -o scf_ex.o
#endif
