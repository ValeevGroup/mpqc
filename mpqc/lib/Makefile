TOPDIR=..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

SUBDIRS = basis elisp perl

include $(SRCDIR)/$(TOPDIR)/lib/GlobalSubDirs

oclean::
	/bin/rm -f lib*.a
	/bin/rm -f lib*.so*
	/bin/rm -f lib*.la
	/bin/rm -rf .libs

clean:: oclean

distclean:: oclean

install:: atominfo.kv magic
	$(INSTALL) $(INSTALLDIROPT) $(scdatadir)
	$(INSTALL) $(INSTALLLIBOPT) $^ $(scdatadir)

