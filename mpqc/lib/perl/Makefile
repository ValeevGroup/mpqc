TOPDIR=../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

install::
	$(INSTALL) $(INSTALLDIROPT) $(libdir)/perl
	$(INSTALL) $(INSTALLLIBOPT) $(wildcard $(SRCDIR)/*.pm) $(libdir)/perl
