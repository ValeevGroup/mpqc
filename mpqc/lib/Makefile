TOPDIR=..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

SUBDIRS = basis elisp perl

include $(SRCDIR)/$(TOPDIR)/lib/GlobalSubDirs

distclean::
	/bin/rm -f lib*.a
	/bin/rm -f lib*.so*

install:: atominfo.kv magic
	$(INSTALL) $(INSTALLDIROPT) $(libdir)
	$(INSTALL) $(INSTALLLIBOPT) $^ $(libdir)

