TOPDIR=.
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

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
