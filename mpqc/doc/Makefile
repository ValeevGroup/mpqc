TOPDIR=..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

default: info

info:
	makeinfo SC.texi

html:
	texi2html -verbose -deepsplit -menu SC.texi

clean distclean::
	/bin/rm -f *.info* *.html

