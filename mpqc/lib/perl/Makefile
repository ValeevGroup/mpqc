TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

install::
	$(INSTALL) $(INSTALLDIROPT) $(scdatadir)/perl
	$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/*.pm $(scdatadir)/perl
