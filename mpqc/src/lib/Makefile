TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

SUBDIRS = util math chemistry

include $(SRCDIR)/$(TOPDIR)/lib/GlobalSubDirs

install::
ifneq ($(LIBTOOL),)
	$(LIBTOOL) --finish $(installroot)$(libdir)
endif

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
