TOPDIR=../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

install::
	$(INSTALL) $(INSTALLDIROPT) $(libdir)/basis
	$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/*.kv $(libdir)/basis
