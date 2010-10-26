TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
include $(TOPDIR)/lib/Makedirlist

SUBDIRS = molecule solvent qc
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_CCA),yes)
SUBDIRS := $(SUBDIRS) cca
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalSubDirs

