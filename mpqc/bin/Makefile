TOPDIR=..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

install::
	$(INSTALL) $(INSTALLDIROPT) $(bindir)
	$(INSTALL) $(INSTALLSCRIPTOPT) sc-config $(bindir)
