TOPDIR=../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

FILES = $(shell echo `find repo`) 

default::

install:: mpqc.cca
	$(INSTALL) $(INSTALLDIROPT) $(installroot)$(libdir)/cca
	$(INSTALL) $(INSTALLIBOPT) mpqc.cca $(installroot)$(libdir)/cca

install_devel:: $(SRCDIR)/sidl/chemistry-mpqc.sidl
	$(INSTALL) $(INSTALLDIROPT) $(installroot)$(libdir)/cca
	$(INSTALL) $(INSTALLDIROPT) $(installroot)$(libdir)/cca/sidl
	$(INSTALL) $(INSTALLLIBOPT) $(SRCDIR)/sidl/chemistry-mpqc.sidl \
          $(installroot)$(libdir)/cca/sidl
	for i in $(FILES); do \
         if test -d $${i}; then \
          $(INSTALL) $(INSTALLDIROPT) $${i} $(installroot)$(libdir)/cca; \
         else \
          $(INSTALL) -D $(INSTALLLIBOPT) $${i} \
            $(installroot)$(libdir)/cca/$${i}; \
         fi; \
        done

clean:: oclean
	if test -d repo; then \
	  rm -rf repo; \
        fi

oclean::
	/bin/rm -f lib*.a
	/bin/rm -f lib*.so*
	/bin/rm -f lib*.la
	/bin/rm -rf .libs

dclean::

distclean:
	/bin/rm -f components.cca

