TOPDIR=..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

DOXYGEN=doxygen

.PHONY: all
all:
	$(DOXYGEN) doxygen.cfg

.PHONY: veryclean
veryclean:
	/bin/rm -rf html latex man
	/bin/rm -rf *~

.PHONY: clean
clean:
	/bin/rm -rf *~

install:
	$(INSTALL) $(INSTALLDIROPT) $(prefix)
	/bin/cp -r html $(prefix)

install_devel:
	$(INSTALL) $(INSTALLDIROPT) $(prefix)
	/bin/cp -r man $(prefix)
	$(INSTALL) $(INSTALLDIROPT) $(prefix)/examples
	$(INSTALL) $(INSTALLDIROPT) $(prefix)/examples/mp2
	$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/devsamp/mp2.cc \
		$(prefix)/examples/mp2
	$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/devsamp/mp2.in \
		$(prefix)/examples/mp2
	sed "s+/usr/local/mpqc/current+$(prefix)+" < \
		$(SRCDIR)/devsamp/Makefile \
		> $(prefix)/examples/mp2/Makefile
