TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

SUBDIRS = util math chemistry

include $(SRCDIR)/$(TOPDIR)/lib/GlobalSubDirs

install::
	$(LIBTOOL) --finish $(installroot)$(libdir)

install_inc::
	$(INSTALL) $(INSTALLDIROPT) $(installroot)$(includedir)
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/*.h $(installroot)$(includedir)
	-$(INSTALL) $(INSTALLLIBOPT) *.h $(installroot)$(includedir)

install_devel::
	$(INSTALL) $(INSTALLDIROPT) $(installroot)$(includedir)
	-$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/*.h $(installroot)$(includedir)
	-$(INSTALL) $(INSTALLLIBOPT) *.h $(installroot)$(includedir)

distclean::
	/bin/rm -f scconfig.h
	/bin/rm -f scdirlist.h
