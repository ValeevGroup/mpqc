TOPDIR=../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

SUBDIRS = util math chemistry

include $(SRCDIR)/$(TOPDIR)/lib/GlobalSubDirs

install:: install_inc

install_inc::
	$(INSTALL) $(INSTALLDIROPT) $(includedir)
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/*.h $(includedir)
	-$(INSTALL) $(INSTALLLIBOPT) *.h $(includedir)
