TOPDIR=..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

install_devel::
	$(INSTALL) $(INSTALLDIROPT) $(bindir)
	$(INSTALL) $(INSTALLSCRIPTOPT) sc-config $(bindir)
	$(INSTALL) $(INSTALLSCRIPTOPT) mkf77sym.pl $(bindir)/sc-mkf77sym
	-$(INSTALL) $(INSTALLSCRIPTOPT) sc-libtool $(bindir)

distclean::
	-/bin/rm -f sc-config
	-/bin/rm -f mkf77sym.pl
