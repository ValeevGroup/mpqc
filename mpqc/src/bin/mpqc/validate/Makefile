TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

INSUF=in
PROGRAM=mpqc
RUN=run
REF=ref

MPQC=../../mpqc

MAKEIN=$(PERL) -I$(SRCTOPDIR)/lib/perl $(SRCDIR)/makein.pl
CHECKOUT=$(PERL) -I$(SRCTOPDIR)/lib/perl $(SRCDIR)/checkout.pl

H2OMASTER=h2o.qci
H2OINPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(RUN) $(H2OMASTER))
H2OOUTPUTS = $(H2OINPUTS:%.$(INSUF)=%.out)
H2OMP2MASTER=h2omp2.qci
H2OMP2INPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(RUN) $(H2OMP2MASTER))
H2OMP2OUTPUTS = $(H2OMP2INPUTS:%.$(INSUF)=%.out)
H2OFRQMASTER=h2ofrq.qci
H2OFRQINPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(RUN) $(H2OFRQMASTER))
H2OFRQOUTPUTS = $(H2OFRQINPUTS:%.$(INSUF)=%.out)
BASIS1MASTER=basis1.qci
BASIS1INPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(RUN) $(BASIS1MASTER))
BASIS1OUTPUTS = $(BASIS1INPUTS:%.$(INSUF)=%.out)
BASIS2MASTER=basis2.qci
BASIS2INPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(RUN) $(BASIS2MASTER))
BASIS2OUTPUTS = $(BASIS2INPUTS:%.$(INSUF)=%.out)
OPTMASTER=opt.qci
OPTINPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(RUN) $(OPTMASTER))
OPTOUTPUTS = $(OPTINPUTS:%.$(INSUF)=%.out)
OPTTSMASTER=optts.qci
OPTTSINPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(RUN) $(OPTTSMASTER))
OPTTSOUTPUTS = $(OPTTSINPUTS:%.$(INSUF)=%.out)

ALLINPUTS = $(H2OINPUTS) $(H2OMP2INPUTS) $(H2OFRQINPUTS) \
            $(BASIS1INPUTS) $(BASIS2INPUTS) \
            $(OPTINPUTS) $(OPTTSINPUTS)
ALLINPUTS = $(H2OINPUTS) $(H2OMP2INPUTS) $(H2OFRQINPUTS) \
            $(BASIS1INPUTS) $(BASIS2INPUTS) \
            $(OPTINPUTS) $(OPTTSINPUTS)
ALLOUTPUTS = $(ALLINPUTS:%.$(INSUF)=%.out)

default::
	@echo \'make inputs\' to make a run directory containing all inputs
	@echo \'make checkrun\' to check outputs of the run directory
	@echo \'make check\' to compare the outputs in run with those in ref
	@echo \'make diff\' to use the diff program to compare outputs

h2o $(H2OINPUTS): $(RUN) $(H2OMASTER)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(RUN) $(H2OMASTER)

h2omp2 $(H2OMP2INPUTS): $(RUN) $(H2OMP2MASTER)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(RUN) $(H2OMP2MASTER)

h2ofrq $(H2OFRQINPUTS): $(RUN) $(H2OFRQMASTER)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(RUN) $(H2OFRQMASTER)

basis1 $(BASIS1INPUTS): $(RUN) $(BASIS1MASTER)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(RUN) $(BASIS1MASTER)

basis2 $(BASIS2INPUTS): $(RUN) $(BASIS2MASTER)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(RUN) $(BASIS2MASTER)

opt $(OPTINPUTS): $(RUN) $(OPTMASTER)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(RUN) $(OPTMASTER)

optts $(OPTTSINPUTS): $(RUN) $(OPTTSMASTER)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(RUN) $(OPTTSMASTER)

$(RUN):
	mkdir -p $(RUN)

inputs: $(ALLINPUTS)

checkrun:
	for file in $(ALLOUTPUTS); \
	  do \
	    ($(CHECKOUT) $${file} ); \
	  done

check:
	for file in $(ALLOUTPUTS:$(RUN)/%=%); \
	  do \
	    ($(CHECKOUT) $(REF)/$${file} $(RUN)/$${file} ); \
	  done

diff:
	-diff -ur $(REF) $(RUN)

$(RUN)/%.out: $(RUN)/%.in
	(cd $(RUN); $(MPQC) -f `basename $<` >& `basename $@`)
