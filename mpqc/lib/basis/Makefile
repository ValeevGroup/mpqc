TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

install::
	$(INSTALL) $(INSTALLDIROPT) $(libdir)/basis
	$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/*.kv $(libdir)/basis

install_devel::
