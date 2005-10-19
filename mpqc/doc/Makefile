TOPDIR=..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

DOXYGEN=doxygen

# Get rid of end-of-comments.  This is done so all pages can be combined
# into a single man page.
MAN_SED_PAGE = (sed s.*/..)

# This makes pages into sections and sections into subsections.  That way
# all the documentation is stuck into a single man page.  The opening
# comment is also removed, because all the pages have to be combined into a
# single comment to suit doxygen.
MAN_SED_SEC = (sed 's.^*/.\<br\>.' | sed 's+\\subsection+\\subsubsection+'| sed 's+\\section+\\subsection+' | sed 's+/\**.*\\page+\\section+')

# This makes the bullets into an 'o', which looks better in the text man pages.
NROFF_SED = sed 's/\\(bu/o/'

.PHONY: all
all:: html
ifeq ($(DOXYGEN_MAN),YES)
all:: man1 man3
endif

.PHONY: man1
man1: mpqc.man.dox
	$(DOXYGEN) doxygen.man1.cfg
	$(NROFF_SED) < man/man1/mpqc.1 > man/man1/mpqc.1.tmp
	/bin/mv man/man1/mpqc.1.tmp man/man1/mpqc.1

.PHONY: man3
man3:
	$(DOXYGEN) doxygen.man3.cfg

.PHONY: mpqc.man.dox
mpqc.man.dox:
	/bin/rm -f $@
	$(MAN_SED_PAGE) < $(SRCDIR)/$(TOPDIR)/src/bin/mpqc/mpqc.dox > $@
	$(MAN_SED_SEC) < $(SRCDIR)/$(TOPDIR)/src/bin/mpqc/mpqcover.dox >> $@
	$(MAN_SED_SEC) < $(SRCDIR)/$(TOPDIR)/src/bin/mpqc/mpqcrunning.dox >> $@
	$(MAN_SED_SEC) < $(SRCDIR)/$(TOPDIR)/src/bin/mpqc/mpqcinp.dox >> $@
	$(MAN_SED_SEC) < $(SRCDIR)/$(TOPDIR)/src/bin/mpqc/mpqcsimp.dox >> $@
	$(MAN_SED_SEC) < $(SRCDIR)/$(TOPDIR)/src/bin/mpqc/mpqcoo.dox >> $@
	$(MAN_SED_SEC) < $(SRCDIR)/$(TOPDIR)/src/bin/mpqc/mpqcval.dox >> $@
	$(MAN_SED_SEC) < $(SRCDIR)/$(TOPDIR)/src/lib/chemistry/qc/psi/mpqcpsi.dox >> $@
	$(MAN_SED_SEC) < $(SRCDIR)/$(TOPDIR)/src/lib/chemistry/cca/mpqccomponents.dox >> $@
	$(MAN_SED_SEC) < $(SRCDIR)/$(TOPDIR)/src/bin/mpqc/mpqclic.dox >> $@
	$(MAN_SED_SEC) < $(SRCDIR)/$(TOPDIR)/src/bin/mpqc/mpqcwar.dox >> $@
	echo "*/" >> $@

.PHONY: html
html:
	$(DOXYGEN) doxygen.cfg

.PHONY: veryclean
veryclean:
	/bin/rm -rf html latex man
	/bin/rm -rf *~

.PHONY: clean
clean:
	/bin/rm -rf *~

install:
	$(INSTALL) $(INSTALLDIROPT) $(installroot)$(prefix)
	/bin/cp -r html $(installroot)$(prefix)

install_devel: install_man install_samples

.PHONY: install_man
install_man:
	$(INSTALL) $(INSTALLDIROPT) $(installroot)$(prefix)
	/bin/cp -r man $(installroot)$(prefix)

.PHONY: install_samples
install_samples:
	$(INSTALL) $(INSTALLDIROPT) $(installroot)$(prefix)
	$(INSTALL) $(INSTALLDIROPT) $(installroot)$(prefix)/examples
	$(INSTALL) $(INSTALLDIROPT) $(installroot)$(prefix)/examples/mp2
	$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/devsamp/mp2.cc \
		$(installroot)$(prefix)/examples/mp2
	$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/devsamp/mp2.in \
		$(installroot)$(prefix)/examples/mp2
	sed "s+/usr/local/mpqc/current+$(prefix)+" < \
		$(SRCDIR)/devsamp/Makefile \
		> $(installroot)$(prefix)/examples/mp2/Makefile
