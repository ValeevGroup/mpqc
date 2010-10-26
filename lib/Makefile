TOPDIR=..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
include $(TOPDIR)/lib/Makedirlist

SUBDIRS = basis elisp perl
ifeq ($(HAVE_SC_LIB_CCA),yes)
SUBDIRS := $(SUBDIRS) cca
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalSubDirs

oclean::
	/bin/rm -f lib*.a
	/bin/rm -f lib*.so*
	/bin/rm -f lib*.la
	/bin/rm -rf .libs

clean:: oclean

distclean:: oclean

install:: atominfo.kv magic
	$(MKINSTALLDIRS) $(installroot)$(scdatadir)
	$(INSTALL) $(INSTALLLIBOPT) $^ $(installroot)$(scdatadir)

