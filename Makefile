TOPDIR=.
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

LOCALMAKEFILE_OPTIONAL = yes

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

ifeq ($(LOCALMAKEFILE_FOUND),yes)

SUBDIRS = lib bin src

include $(SRCDIR)/$(TOPDIR)/lib/GlobalSubDirs

install_devel::
	-$(INSTALL) $(INSTALLSCRIPTOPT) libtool $(installroot)$(bindir)/sc-libtool

clean::
	-rm -f depcheck.cc

distclean:: clean
	-rm -f config.log
	-rm -f config.status
	-rm -f libtool

else

# If the LocalMakefile does not exist we are not in an object directory.
# Here we only wish to do administrative tasks, such as make tags.

default::
	@echo "The LocalMakefile was not found.  First run configure in the directory"
	@echo "you want object code to be placed and then run make there.  However,"
	@echo "if you would like to run some administrative commands in the source"
	@echo "tree, the following targets are available:"
	@echo 
	@echo "  configure     Update aclocal.m4 and run configure."
	@echo "  touch         Make sure the parser files are more recent than"
	@echo "                the sources.  Useful after a CVS checkout."
	@echo "  tags          Build a new TAGS file."
	@echo "  ebrowse       Build a new BROWSE file for emacs ebrowse."
	@echo "  clean         Remove emacs backup files in this and all subdirectories."
	@echo 

endif

check: check0

check0:
	cd src/bin/mpqc/validate; $(MAKE) check0

check1:
	cd src/bin/mpqc/validate; $(MAKE) check1

check2:
	cd src/bin/mpqc/validate; $(MAKE) check2

check_clean:
	cd src/bin/mpqc/validate; $(MAKE) check_clean

.PHONY: configure
configure:
	aclocal -I lib/autoconf
	autoconf
	/bin/rm -rf autom4te.cache

.PHONY: touch
touch:
	touch src/bin/mpqc/scan.cc
	touch src/bin/mpqc/parse.cc
	touch src/bin/mpqc/parse.h
	touch src/lib/util/keyval/ipv2_scan.cc
	touch src/lib/util/keyval/ipv2_parse.cc
	touch src/lib/util/keyval/ipv2_parse.h

.PHONY: tags
tags:
	etags --members `find . -name "*.[hcfCF]"` `find . -name "*.cc"`

.PHONY: ebrowse
ebrowse:
	ebrowse `find . -name "*.[hcC]"` `find . -name "*.cc"`

.PHONY: clean
clean::
	/bin/rm -f `find . -name "*~" -print`
