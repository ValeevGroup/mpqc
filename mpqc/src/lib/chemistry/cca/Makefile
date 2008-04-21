TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

default::

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
include $(TOPDIR)/lib/Makedirlist
include $(SRCDIR)/$(TOPDIR)/lib/cca/CCAMakefile

ifneq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_CINTS),yes)
  EXCLUDE = --exclude='.*CintsEvaluatorFactory.*'
endif
ifneq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_LIBINT2),yes)
  EXCLUDE = --exclude='.*Libint2EvaluatorFactory.*'
endif

SUBDIRS = server int molecule

ifeq ($(ENABLESHARED),yes)
  SUBDIRS += module
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalSubDirs

