TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

SUBDIRS = util math chemistry

include $(SRCDIR)/$(TOPDIR)/lib/GlobalSubDirs

install::

install_inc::
	$(INSTALL) $(INSTALLDIROPT) $(includedir)
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/*.h $(includedir)
	-$(INSTALL) $(INSTALLLIBOPT) *.h $(includedir)

install_devel::
	$(INSTALL) $(INSTALLDIROPT) $(includedir)
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/*.h $(includedir)
	-$(INSTALL) $(INSTALLLIBOPT) *.h $(includedir)
