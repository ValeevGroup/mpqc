TOPDIR=..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

T2H = $(SRCDIR)/$(TOPDIR)/bin/texi2html

TEXIDIRS = $(SRCDIR)
EXTRATEXINPUT = $(shell /bin/sh -f $(SRCDIR)/addcolon.sh $(TEXIDIRS))
INCLUDEDIRS = $(shell /bin/sh -f $(SRCDIR)/addi.sh $(TEXIDIRS))

default: info

.PHONY:	info html dvisee pssee SC.dvi

info: SC.texi
	makeinfo $(INCLUDEDIRS) $^

html: SC.texi
	TEXINPUTS=$(EXTRATEXINPUT):$$TEXINPUTS & perl $(T2H) -verbose -deepsplit -menu $^

SC.dvi: SC.texi
	TEXINPUTS=$(EXTRATEXINPUT):$$TEXINPUTS & texi2dvi $^

SC.ps: SC.dvi
	dvips -o SC.ps $^

dvisee: SC.dvi
	xdvi $^

pssee: SC.ps
	ghostview $^

echo:
	echo $(EXTRATEXINPUT):$(TEXINPUTS)
	echo $(INCLUDEDIRS)

clean distclean::
	/bin/rm -f *.info* *.html
	/bin/rm -f SC.{aux,cl,cls,cp,cps,dvi,fn,ky,log,mc,mcs}
	/bin/rm -f SC.{pg,toc,tp,vr,ps}
	/bin/rm -f SC.*.bak


