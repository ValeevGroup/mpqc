TOPDIR=..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

T2H = $(SRCDIR)/$(TOPDIR)

default: info

.PHONY:	info html

info: SC.texi
	makeinfo $^

html: SC.texi
	perl $(T2H)/bin/texi2html -verbose -deepsplit -menu $^

clean distclean::
	/bin/rm -f *.info* *.html

