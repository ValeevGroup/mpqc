TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TEXINPUTS=$(SRCDIR)/$(TOPDIR)/src/lib:$(SRCDIR):$(SRCDIR)/$(TOPDIR)/lib/texinputs:$(TOPDIR)/src/lib
LATEX=latex
XDVI=xdvi
LATEX2HTML=latex2html
L2HFLAGS=-dir /net/munin/usr/people/cljanss/public_html/mpqc/user \
         -up_url "http://midway.ca.sandia.gov/%7Ecljanss/mpqc.html" \
         -up_title "The MPQC Home Page" -html_version 3.1 \
         -split 4 -link 4 \
         -index_in_navigation -contents_in_navigation -bottom_navigation

default:: mpqcic.see

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

mpqcic.see: mpqcic.dvi
	$(XDVI) $^

mpqcic.dvi: mpqcic.tex
	TEXINPUTS=$$TEXINPUTS:$(TEXINPUTS) $(LATEX) $^
	makeindex mpqcic.idx
	TEXINPUTS=$$TEXINPUTS:$(TEXINPUTS) $(LATEX) $^
			      
mpqcic.html: mpqcic.tex	      
	TEXINPUTS=$$TEXINPUTS:$(TEXINPUTS) $(LATEX2HTML) $(L2HFLAGS) $^