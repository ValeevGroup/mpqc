TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

# the $(shell pwd) is needed for latex2html
EXTRA_TEXINPUTS=$(SRCTOPDIR)/src/lib:$(SRCDIR):$(SRCTOPDIR)/lib/texinputs:$(TOPDIR)/src/lib:$(shell pwd)
LATEX=latex
XDVI=xdvi
LATEX2HTML=latex2html
L2HFLAGS=-dir mpqcuser \
         -up_url ".." \
         -up_title "The MPQC Home Page" -html_version 3.2 \
         -split 4 -link 4 -local_icons -no_math -no_reuse \
         -index_in_navigation -contents_in_navigation -bottom_navigation

SAMPLEINPUTS = scf.in.tex scfopt.in.tex mp2.in.tex scffreq.in.tex \
               hbondopt.in.tex mancoor.in.tex fixed.in.tex ts.in.tex \
               scfoptguesshess.in.tex tsguesshess.in.tex newton.in.tex

DOCFILES = classhier.tex compiling.tex input.tex mpqc.tex running.tex \
           version.tex

default:: mpqc.see

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

mpqc.see: mpqc.dvi
	$(XDVI) $^

mpqc.ps: mpqc.dvi
	dvips -o $@ $^

mpqc.dvi: mpqc.tex $(SAMPLEINPUTS) $(DOCFILES)
	@echo env TEXINPUTS=$$TEXINPUTS
	@echo TEXINPUTS=$$TEXINPUTS:$(EXTRA_TEXINPUTS)
	TEXINPUTS=$$TEXINPUTS:$(EXTRA_TEXINPUTS) $(LATEX) $<
	makeindex mpqc.idx
	TEXINPUTS=$$TEXINPUTS:$(EXTRA_TEXINPUTS) $(LATEX) $<
			      
mpqc.html: mpqc.tex $(SAMPLEINPUTS) $(DOCFILES)
	TEXINPUTS=$$TEXINPUTS:$(EXTRA_TEXINPUTS) $(LATEX2HTML) $(L2HFLAGS) $<

targetclean::
	/bin/rm -rf mpqc.dvi mpqc.html

%.in.tex: $(SRCDIR)/../sample/%.in
	$(PERL) $(SRCDIR)/mpqctotex.pl $^ > $@

version.tex: $(SRCDIR)/../version.h
	sed -e "s/.define MPQC_VERSION //" -e "s/\"//g" < $< > $@
