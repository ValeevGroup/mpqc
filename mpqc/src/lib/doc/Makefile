TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

# The PWD is needed for a latex2html bug
EXTRA_TEXINPUTS=$(SRCTOPDIR)/src/lib:$(SRCDIR):$(SRCTOPDIR)/lib/texinputs:$(TOPDIR)/src/lib:$(shell pwd):$(PWD)/$(TOPDIR)/src/lib
#EXTRA_TEXINPUTS=$(SRCDIR)/$(TOPDIR)/src/lib:$(SRCDIR):$(SRCDIR)/$(TOPDIR)/lib/texinputs:$(PWD)/$(TOPDIR)/src/lib
LATEX2HTMLSTYLES=$(SRCDIR)/$(TOPDIR)/lib/texinputs
LATEX=latex
XDVI=xdvi
LATEX2HTML=latex2html
L2HFLAGS=-dir scprog \
         -up_url ".." \
         -up_title "The SC Toolkit Home Page" -html_version 3.2 \
         -split 4 -link 4 -local_icons -no_reuse \
         -index_in_navigation -contents_in_navigation -bottom_navigation

L2HFLAGSTEST= -dir /net/munin/usr/people/cljanss/public_html/mpqc/test \
              -split 4 -link 4 \
              -index_in_navigation -contents_in_navigation -bottom_navigation

default:: prog.see

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

%.see: %.dvi
	$(XDVI) $^

%.dvi: %.tex
	@echo env TEXINPUTS=$$TEXINPUTS
	@echo TEXINPUTS=$$TEXINPUTS:$(EXTRA_TEXINPUTS)
	TEXINPUTS=$$TEXINPUTS:$(EXTRA_TEXINPUTS) $(LATEX) $<
	makeindex $(*F).idx
	TEXINPUTS=$$TEXINPUTS:$(EXTRA_TEXINPUTS) $(LATEX) $<

%.html: %.tex #$(SRCDIR)/scdoc.perl
	TEXINPUTS=$$TEXINPUTS:$(EXTRA_TEXINPUTS) $(LATEX2HTML) $(L2HFLAGS) $<

# this is needed for latex2html
#$(SRCDIR)/scdoc.perl:
#	cd $(SRCDIR) && ln -s ../../../lib/texinputs/scdoc.perl scdoc.perl
