TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

install::
	$(INSTALL) $(INSTALLDIROPT) $(scdatadir)/basis
	$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/*.kv $(scdatadir)/basis

install_devel::
