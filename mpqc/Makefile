TOPDIR=.
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

SUBDIRS = lib bin src

include $(SRCDIR)/$(TOPDIR)/lib/GlobalSubDirs

miscclean::
	-rm -f config.log
