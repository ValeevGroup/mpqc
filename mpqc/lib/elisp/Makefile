TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

install::
	$(INSTALL) $(INSTALLDIROPT) $(scdatadir)/elisp
	$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/*.el $(scdatadir)/elisp
