TOPDIR=..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

C2T = $(SRCDIR)/$(TOPDIR)/bin/class2texi
T2H = $(SRCDIR)/$(TOPDIR)/bin/texi2html
RMDOTDOT = $(SRCDIR)/$(TOPDIR)/bin/rmdotdot

SRCLIBDIR := $(shell $(RMDOTDOT) $(SRCDIR)/$(TOPDIR)/src.lib)

# The are the hand written texi files.
DOC =  SC.texi class.texi idxcon.texi keyval.texi scmat.texi \
       idxcls.texi idxmac.texi ref.texi state.texi

# These are the files generated from keyval.h
KEYVALDOC = KeyVal.texi KeyVal.pub.texi KeyVal.pro.texi

# These are all of texi files used to generate the documentation.
# The top level document, SC.texi, must be first in the list.
ALLDOC = $(DOC) $(KEYVALDOC)

TEXIDIRS = $(SRCDIR)
EXTRATEXINPUT = $(shell /bin/sh -f $(SRCDIR)/addcolon.sh $(TEXIDIRS))
INCLUDEDIRS = $(shell /bin/sh -f $(SRCDIR)/addi.sh $(TEXIDIRS))

default: info

.PHONY:	info html dvisee pssee SC.dvi

info: $(ALLDOC)
	makeinfo $(INCLUDEDIRS) $<

html: $(ALLDOC)
	TEXINPUTS=$(EXTRATEXINPUT):$$TEXINPUTS & perl $(T2H) -verbose -deepsplit -menu $<

SC.dvi: $(ALLDOC)
	TEXINPUTS=$(EXTRATEXINPUT):$$TEXINPUTS & texi2dvi $<

SC.ps: SC.dvi
	dvips -o SC.ps $<

dvisee: SC.dvi
	xdvi $<

pssee: SC.ps
	ghostview $<

echo:
	echo $(EXTRATEXINPUT):$(TEXINPUTS)
	echo $(INCLUDEDIRS)

clean distclean::
	/bin/rm -f *.info* *.html
	/bin/rm -f SC.{aux,cl,cls,cp,cps,dvi,fn,ky,log,mc,mcs}
	/bin/rm -f SC.{pg,toc,tp,vr,ps}
	/bin/rm -f SC.*.bak

$(KEYVALDOC): $(SRCLIBDIR)/util/keyval/keyval.h
	$(C2T) $<
