TOPDIR=..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

install_devel::
	$(MKINSTALLDIRS) $(installroot)$(bindir)
	$(INSTALL) $(INSTALLSCRIPTOPT) sc-config $(installroot)$(bindir)
	$(INSTALL) $(INSTALLSCRIPTOPT) mkf77sym.pl $(installroot)$(bindir)/sc-mkf77sym

distclean::
	-/bin/rm -f sc-config
	-/bin/rm -f mkf77sym.pl
