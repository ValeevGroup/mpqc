TOPDIR=../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

install::
	$(INSTALL) $(INSTALLDIROPT) $(libdir)/elisp
	$(INSTALL) $(INSTALLLIBOPT) $(wildcard $(SRCDIR)/*.el) $(libdir)/elisp
