TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

# The PWD is needed for a latex2html bug
TEXINPUTS=$(SRCDIR)/$(TOPDIR)/src/lib:$(SRCDIR):$(SRCDIR)/$(TOPDIR)/lib/texinputs:$(PWD)/$(TOPDIR)/src/lib
LATEX2HTMLSTYLES=$(SRCDIR)/$(TOPDIR)/lib/texinputs
LATEX=latex
XDVI=xdvi
LATEX2HTML=latex2html
L2HFLAGS=-dir /net/munin/usr/people/cljanss/public_html/mpqc/prog.new \
         -up_url "http://midway.ca.sandia.gov/%7Ecljanss/mpqc.html" \
         -up_title "The MPQC Home Page" -html_version 3.1 \
         -split 4 -link 4 \
         -index_in_navigation -contents_in_navigation -bottom_navigation

L2HFLAGSTEST= -dir /net/munin/usr/people/cljanss/public_html/mpqc/test \
              -split 4 -link 4 \
              -index_in_navigation -contents_in_navigation -bottom_navigation

default:: prog.see

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

%.see: %.dvi
	$(XDVI) $^

%.dvi: %.tex
	TEXINPUTS=$$TEXINPUTS:$(TEXINPUTS) $(LATEX) $^
	makeindex $(*F).idx
	TEXINPUTS=$$TEXINPUTS:$(TEXINPUTS) $(LATEX) $^

%.html: %.tex $(SRCDIR)/scdoc.perl
#	export LATEX2HTMLSTYLES; export TEXINPUTS; \
#	LATEX2HTMLSTYLES=$(LATEX2HTMLSTYLES) && TEXINPUTS=$$TEXINPUTS:$(TEXINPUTS) && $(LATEX2HTML) $(L2HFLAGS) $<
	export TEXINPUTS; \
	TEXINPUTS=$$TEXINPUTS:$(TEXINPUTS) && $(LATEX2HTML) $(L2HFLAGS) $<

# this is needed for latex2html
$(SRCDIR)/scdoc.perl:
	cd $(SRCDIR) && ln -s ../../../lib/texinputs/scdoc.perl scdoc.perl
