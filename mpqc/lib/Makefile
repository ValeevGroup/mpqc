TOPDIR=..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

SUBDIRS = basis elisp perl texinputs

include $(SRCDIR)/$(TOPDIR)/lib/GlobalSubDirs

distclean::
	/bin/rm -f lib*.a
	/bin/rm -f lib*.so*

install:: atominfo.kv magic
	$(INSTALL) $(INSTALLDIROPT) $(libdir)
	$(INSTALL) $(INSTALLLIBOPT) $^ $(libdir)

