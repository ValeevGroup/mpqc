TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
include $(TOPDIR)/lib/Makedirlist

INSUF=in
PROGRAM=mpqc
RUN=run
INP=input
REF=$(SRCDIR)/ref
MPQC=../../mpqc

INPUTGENDEPS=$(SRCTOPDIR)/lib/perl/QCParse.pm \
             $(SRCTOPDIR)/lib/perl/QCResult.pm \
             $(SRCTOPDIR)/lib/perl/Molecule.pm \
             $(SRCDIR)/makein.pl

MAKEIN=$(PERL) -I$(SRCTOPDIR)/lib/perl $(SRCDIR)/makein.pl
CHECKOUT=$(PERL) -I$(SRCTOPDIR)/lib/perl $(SRCDIR)/checkout.pl

########################################################################
# need to determine nthreadperproc from MPQCRUN_ARGS, since this determines
# whether or not certain tests

ifneq ($(filter --nthreadperproc,$(MPQCRUN_ARGS)),)
$(error must use --nthreadperproc=xxx syntax for MPQCRUN_ARGS)
endif
MPQCRUN_NTHREADPERPROC:=$(filter --nthreadperproc=%, $(MPQCRUN_ARGS))
MPQCRUN_NTHREADPERPROC:=$(patsubst --nthreadperproc=%, %, $(MPQCRUN_NTHREADPERPROC))
MPQCRUN_NTHREADPERPROC:=$(strip $(MPQCRUN_NTHREADPERPROC))
ifeq ($(MPQCRUN_NTHREADPERPROC),)
  MPQCRUN_NTHREADPERPROC=1
endif

ifneq ($(filter --integral,$(MPQCRUN_ARGS)),)
$(error must use --integral=xxx syntax for MPQCRUN_ARGS)
endif
MPQCRUN_INTEGRAL:=$(filter --integral=%, $(MPQCRUN_ARGS))
MPQCRUN_INTEGRAL:=$(patsubst --integral=%, %, $(MPQCRUN_INTEGRAL))
MPQCRUN_INTEGRAL:=$(strip $(MPQCRUN_INTEGRAL))
ifeq ($(MPQCRUN_INTEGRAL),)
  MPQCRUN_INTEGRAL=intv3
endif

# Also need to see if --mpqc was specified.  If not, then add
# --mpqc ../mpqc.
ifeq ($(filter --mpqc, $(MPQCRUN_ARGS)),)
MPQCRUN_ARGS:=--mpqc ../mpqc $(MPQCRUN_ARGS)
endif

ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_MBPTR12),yes)
ifeq ($(MPQCRUN_INTEGRAL),cints)
HAVE_MBPTR12=yes
endif
endif

########################################################################

H2OMASTER=h2o.qci
H2OINPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(INP) $(H2OMASTER))
H2OOUTPUTS = $(H2OINPUTS:%.$(INSUF)=%.out)

H2OMP2MASTER=h2omp2.qci
H2OMP2INPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(INP) $(H2OMP2MASTER))
H2OMP2OUTPUTS = $(H2OMP2INPUTS:%.$(INSUF)=%.out)

MP2R12MASTER=mp2r12.qci
ifeq ($(HAVE_MBPTR12),yes)
  MP2R12INPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(INP) $(MP2R12MASTER))
else
  MP2R12INPUTS=
endif
MP2R12OUTPUTS = $(MP2R12INPUTS:%.$(INSUF)=%.out)

H2OFRQMASTER=h2ofrq.qci
H2OFRQINPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(INP) $(H2OFRQMASTER))
H2OFRQOUTPUTS = $(H2OFRQINPUTS:%.$(INSUF)=%.out)

ORTHOGMASTER=orthog.qci
ORTHOGINPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(INP) $(ORTHOGMASTER))
ORTHOGOUTPUTS = $(ORTHOGINPUTS:%.$(INSUF)=%.out)

BASIS1MASTER=basis1.qci
BASIS1INPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(INP) $(BASIS1MASTER))
BASIS1OUTPUTS = $(BASIS1INPUTS:%.$(INSUF)=%.out)

DFTMASTER=dft.qci
DFTINPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(INP) $(DFTMASTER))
DFTOUTPUTS = $(DFTINPUTS:%.$(INSUF)=%.out)

BASIS2MASTER=basis2.qci
BASIS2INPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(INP) $(BASIS2MASTER))
BASIS2OUTPUTS = $(BASIS2INPUTS:%.$(INSUF)=%.out)

OPTMASTER=opt.qci
OPTINPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(INP) $(OPTMASTER))
OPTOUTPUTS = $(OPTINPUTS:%.$(INSUF)=%.out)
OPTTSMASTER=optts.qci
OPTTSINPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(INP) $(OPTTSMASTER))
OPTTSOUTPUTS = $(OPTTSINPUTS:%.$(INSUF)=%.out)

SYMM1MASTER=symm1.qci
SYMM1INPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(INP) $(SYMM1MASTER))
SYMM1OUTPUTS = $(SYMM1INPUTS:%.$(INSUF)=%.out)

SYMM2MASTER=symm2.qci
SYMM2INPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(INP) $(SYMM2MASTER))
SYMM2OUTPUTS = $(SYMM2INPUTS:%.$(INSUF)=%.out)

SYMM3MASTER=symm3.qci
SYMM3INPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(INP) $(SYMM3MASTER))
SYMM3OUTPUTS = $(SYMM3INPUTS:%.$(INSUF)=%.out)

CLSCFMASTER=clscf.qci
CLSCFINPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(INP) $(CLSCFMASTER))
CLSCFOUTPUTS = $(CLSCFINPUTS:%.$(INSUF)=%.out)

USCFMASTER=uscf.qci
USCFINPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(INP) $(USCFMASTER))
USCFOUTPUTS = $(USCFINPUTS:%.$(INSUF)=%.out)

HSOSSCFMASTER=hsosscf.qci
HSOSSCFINPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(INP) $(HSOSSCFMASTER))
HSOSSCFOUTPUTS = $(HSOSSCFINPUTS:%.$(INSUF)=%.out)

CCAINTMASTER=ccaint.qci
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_INTCCA),yes)
 CCAINTINPUTS=$(shell $(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -e -d $(INP) $(CCAINTMASTER))
endif
CCAINTOUTPUTS = $(CCAINTINPUTS:%.$(INSUF)=%.out)

CKPTINPUTS = \
 ckpt_0clscf.in     ckpt_0qnewtopt.in  ckpt_1zapt2.in     ckpt_mp2.in \
 ckpt_0efcopt.in    ckpt_1clscf.in     ckpt_1qnewtopt.in  ckpt_2qnewtopt.in \
 ckpt_0hsosscf.in   ckpt_1efcopt.in    ckpt_clscf.in      ckpt_efcopt.in \
 ckpt_0mp2.in       ckpt_1hsosscf.in   ckpt_2efcopt.in    ckpt_qnewtopt.in \
 ckpt_0zapt2.in     ckpt_1mp2.in       ckpt_hsosscf.in \
 ckpt_0clksb3lyp.in     ckpt_1clksb3lyp.in     ckpt_clksb3lyp.in \
 ckpt_0clksbp86.in      ckpt_1clksbp86.in      ckpt_clksbp86.in \
 ckpt_0clkshfg96.in     ckpt_1clkshfg96.in     ckpt_clkshfg96.in \
 ckpt_0clksmpwpw91.in   ckpt_1clksmpwpw91.in   ckpt_clksmpwpw91.in \
 ckpt_0clkspbe.in       ckpt_1clkspbe.in       ckpt_clkspbe.in \
 ckpt_0clkspw91.in      ckpt_1clkspw91.in      ckpt_clkspw91.in \
 ckpt_0clksspz81.in     ckpt_1clksspz81.in     ckpt_clksspz81.in \
 ckpt_0clkssvwn1.in     ckpt_1clkssvwn1.in     ckpt_clkssvwn1.in \
 ckpt_0clkssvwn1rpa.in  ckpt_1clkssvwn1rpa.in  ckpt_clkssvwn1rpa.in \
 ckpt_0clkssvwn2.in     ckpt_1clkssvwn2.in     ckpt_clkssvwn2.in \
 ckpt_0clkssvwn3.in     ckpt_1clkssvwn3.in     ckpt_clkssvwn3.in \
 ckpt_0clkssvwn4.in     ckpt_1clkssvwn4.in     ckpt_clkssvwn4.in \
 ckpt_0clkssvwn5.in     ckpt_1clkssvwn5.in     ckpt_clkssvwn5.in \
 ckpt_0clksxalpha.in    ckpt_1clksxalpha.in    ckpt_clksxalpha.in \
 ckpt_0hsosksxalpha.in  ckpt_1hsosksxalpha.in  ckpt_hsosksxalpha.in \
 ckpt_0uksxalpha.in     ckpt_1uksxalpha.in     ckpt_uksxalpha.in

CKPTINPUTS := $(CKPTINPUTS:%=$(INP)/%)
CKPTOUTPUTS = $(CKPTINPUTS:%.$(INSUF)=%.out)

MBPTINPUTS = \
 mbpt_mp2mem.in     mbpt_opt12v1.in    mbpt_opt22v2.in    mbpt_zapt2v2lb.in \
 mbpt_mp2v1.in      mbpt_opt12v2.in    mbpt_opt22v2lb.in \
 mbpt_mp2v2.in      mbpt_opt12v2lb.in  mbpt_zapt2v1.in \
 mbpt_mp2v2lb.in    mbpt_opt22v1.in    mbpt_zapt2v2.in \
 mbpt_mp2mem_mp.in  mbpt_mp2v2_mp.in   mbpt_mp2mem_c1.in  mbpt_mp2mem_auto.in

ifeq ($(HAVE_MBPTR12),yes)
 MBPTINPUTS += \
    mbpt_mp2r12_ne2.in            mbpt_mp2r12_ne2_posix.in \
    mbpt_mp2r12_ne2_tz.in \
    mbpt_mp2r12_c6h6_multipass.in
 ifeq ($(MPQCRUN_NTHREADPERPROC),1)
  MBPTINPUTS += mbpt_mp2r12_ne2_multipass.in
 endif
endif

MBPTINPUTS := $(MBPTINPUTS:%=$(INP)/%)
MBPTOUTPUTS = $(MBPTINPUTS:%.$(INSUF)=%.out)

METHODSINPUTS = \
clhf.in        hsosks_b3lyp.in       tchf.in \
clks_b3lyp.in  hsosks_b88.in         uhf.in \
clks_b88.in    hsosks_blyp.in        uks_b3lyp.in \
clks_blyp.in   hsosks_lsdax.in       uks_b88.in \
clks_lsdax.in  hsosks_xa.in          uks_blyp.in \
clks_xa.in     uks_lsdax.in          uks_xa.in \
hsoshf.in      osshf.in

ifeq ($(HAVE_MBPTR12),yes)
 METHODSINPUTS += \
    mp2r12ap_+gbc.in mp2r12ap_abs+.in \
    mp2r12ap_abs.in  mp2r12ap_cabs+.in \
    mp2r12ap_cabs.in mp2r12ap_nogebc.in \
    mp2r12ap_+ebc.in mp2r12ap_+gbc+ebc.in
endif

METHODSINPUTS := $(METHODSINPUTS:%=$(INP)/methods_%)

INPUTINPUTS = \
rksch2.in uksch2.in uhfch2opt.in rhfch2opt.in hfch2opt.in \
hfh2ofreq.in ksh2oco.in ksh2o.in hfh2oopt.in \
mp2h2o.in zapt2ch2.in
ifeq ($(HAVE_MBPTR12),yes)
  INPUTINPUTS += mp2r12ah2o.in mp2r12aph2o.in
endif
INPUTINPUTS := $(INPUTINPUTS:%=$(INP)/input_%)

##############################################################################

# cints will not run all test cases.  This filter can be used to select
# just cases that cints can run
CINTSFILTER = $(INP)/mbpt_mp2r12% $(INP)/methods_mp2r12% $(INP)/mp2r12%

ALLINPUTS = $(METHODSINPUTS) $(H2OINPUTS) $(H2OMP2INPUTS) $(H2OFRQINPUTS) \
	    $(MP2R12INPUTS) \
            $(ORTHOGINPUTS) $(BASIS1INPUTS) $(BASIS2INPUTS) \
            $(DFTINPUTS) \
            $(OPTINPUTS) $(OPTTSINPUTS) \
            $(SYMM1INPUTS) $(SYMM2INPUTS) $(SYMM3INPUTS) \
            $(CKPTINPUTS) $(MBPTINPUTS) \
            $(CLSCFINPUTS) $(USCFINPUTS) $(HSOSSCFINPUTS) \
            $(INPUTINPUTS) $(CCAINTINPUTS)
ALLOUTPUTS = $(ALLINPUTS:%.$(INSUF)=%.out)

CHECK0INPUTS := $(H2OINPUTS) $(H2OFRQINPUTS) $(H2OMP2INPUTS) \
                $(METHODSINPUTS) $(INPUTINPUTS)
CHECK0BFILTER := %ccpvdz %ccpvtz %ccpvqz %ccpv5z
CHECK0BFILTER += %ccpcvdz %ccpcvtz %ccpcvqz %ccpcv5z
CHECK0BFILTER += %pc2 %pc3 %pc4
CHECK0BFILTER += %6311gss
CHECK0FILTER := $(CHECK0BFILTER:%=%d2h) \
                $(CHECK0BFILTER:%=%c2v) \
                $(CHECK0BFILTER:%=%c2) \
                $(CHECK0BFILTER:%=%ci) \
                $(CHECK0BFILTER:%=%auto) \
                $(CHECK0BFILTER:%=%cs) \
                $(CHECK0BFILTER:%=%c1)
CHECK0FILTER := $(CHECK0FILTER:%=%.in) \
                $(CHECK0FILTER:%=%frq.in) \
                $(CHECK0FILTER:%= %opt.in) \
                $(CHECK0FILTER:%= %optfrq.in)
CHECK0FILTER += $(INP)/basis2_% $(INP)/dft_% $(INP)/symm1_cub%
CHECK0FILTER += $(INP)/mbpt_mp2r12_c6h6_multipass.in
CHECK0INPUTS := $(filter-out $(CHECK0FILTER),$(CHECK0INPUTS))
ifeq ($(MPQCRUN_INTEGRAL),cints)
  CHECK0INPUTS := $(filter $(CINTSFILTER),$(CHECK0INPUTS))
endif
CHECK0OUTPUTS := $(CHECK0INPUTS:%.$(INSUF)=%.out)
#printcheck0:
#	@echo CHECK0BFILTER = $(CHECK0BFILTER)
#	@echo CHECK0FILTER = $(CHECK0FILTER)
#	@ls -l $(CHECK0INPUTS)

CHECK1INPUTS := $(ALLINPUTS)
CHECK1BFILTER := %ccpvdz %ccpvtz %ccpvqz %ccpv5z
CHECK1BFILTER += %ccpcvdz %ccpcvtz %ccpcvqz %ccpcv5z
CHECK1BFILTER += %pc2 %pc3 %pc4
CHECK1BFILTER += %6311gss
CHECK1FILTER := $(CHECK1BFILTER:%=%d2h) \
                $(CHECK1BFILTER:%=%c2v) \
                $(CHECK1BFILTER:%=%c2) \
                $(CHECK1BFILTER:%=%ci) \
                $(CHECK1BFILTER:%=%auto) \
                $(CHECK1BFILTER:%=%cs) \
                $(CHECK1BFILTER:%=%c1)
CHECK1FILTER := $(CHECK1FILTER:%=%.in) \
                $(CHECK1FILTER:%=%frq.in) \
                $(CHECK1FILTER:%= %opt.in) \
                $(CHECK1FILTER:%= %optfrq.in)
CHECK1FILTER += $(INP)/basis2_% $(INP)/dft_% $(INP)/symm1_cub%
CHECK1FILTER += $(INP)/mbpt_mp2r12_c6h6_multipass.in
CHECK1INPUTS := $(filter-out $(CHECK1FILTER),$(CHECK1INPUTS))
ifeq ($(MPQCRUN_INTEGRAL),cints)
  CHECK1INPUTS := $(filter $(CINTSFILTER),$(CHECK1INPUTS))
endif
CHECK1OUTPUTS := $(CHECK1INPUTS:%.$(INSUF)=%.out)
#printcheck1:
#	@echo CHECK1BFILTER = $(CHECK1BFILTER)
#	@echo CHECK1FILTER = $(CHECK1FILTER)
#	@ls -l $(CHECK1INPUTS)

CHECK2INPUTS := $(ALLINPUTS)
ifeq ($(MPQCRUN_INTEGRAL),cints)
  CHECK2INPUTS := $(filter $(CINTSFILTER),$(CHECK2INPUTS))
endif
CHECK2OUTPUTS := $(CHECK2INPUTS:%.$(INSUF)=%.out)
#printcheck2:
#	@ls -l $(CHECK2INPUTS)

##############################################################################

.PHONY: default
default::
	@echo \'make check0\' to run and check the small test suite
	@echo \'make check0_run\' to run the small test suite
	@echo \'make check0_chk\' to check the small test suite
	@echo \'make check1\' to run and check the intermediate test suite
	@echo \'make check1_run\' to run the small test suite
	@echo \'make check1_chk\' to check the small test suite
	@echo \'make check2\' to inputs and check the full test suite
	@echo \'make check2_run\' to run the small test suite
	@echo \'make check2_chk\' to check the small test suite
	@echo "  For each of the above, MPQCRUN_ARGS=... can be given to"
	@echo "  to control how the jobs are run.  If MPQCRUN_ARGS is"
	@echo "  given to checkn_run, then the same MPQCRUN_ARGS must be given"
	@echo "  to the checkn_chk."
	@echo \'make inputs\' to make a run directory containing all inputs.
	@echo "               This is not needed to run the checks since the"
	@echo "               inputs from the src directory are used.  It"
	@echo "               is only for maintainer use."
	@echo \'make check_clean\' removes output and scratch files from the run directory
	@echo \'make check_clean_scratch\' removes scratch files from the run directory
	@echo Deprecated make targets:
	@echo \'make checkrun\' to check outputs of the run directory
	@echo \'make check\' to compare the outputs in run with those in ref
	@echo \'make diff\' to use the diff program to compare outputs

.PHONY: check check0 check0_run check0_chk
.PHONY: check1 check1_run check1_chk
.PHONY: check2 check2_run check2_chk
.PHONY: check_clean check_clean_scratch

check_clean: check_clean_scratch
	/bin/rm -f $(RUN)/*.out

check_clean_scratch:
	/bin/rm -f $(RUN)/*.tmp $(RUN)/*.wfn $(RUN)/*.hess

check: check0

check0: $(RUN) check0_run check0_chk

MPQCRUN = ../mpqcrun
MPQCRUN_EXTRA_ARGS = --autoout --simpout \
                     --inputprefix $(SRCDIR)/ref/ \
                     --outputprefix $(RUN)/ \
                     --seq ^ckpt_

check0_run:
	$(PERL) $(MPQCRUN) $(MPQCRUN_EXTRA_ARGS) $(MPQCRUN_ARGS) $(CHECK0INPUTS:$(INP)/%=%)

check0_chk:
	$(CHECKOUT) -r $(SRCDIR)/ref $(CHECK0OUTPUTS:$(INP)/%=$(RUN)/%)

check1: $(RUN) check1_run check1_chk

check1_run:
	$(PERL) $(MPQCRUN) $(MPQCRUN_EXTRA_ARGS) $(MPQCRUN_ARGS) $(CHECK1INPUTS:$(INP)/%=%)

check1_chk:
	$(CHECKOUT) -r $(SRCDIR)/ref $(CHECK1OUTPUTS:$(INP)/%=$(RUN)/%)

check2: (RUN) check2_run check2_chk

check2_run:
	$(PERL) $(MPQCRUN) $(MPQCRUN_EXTRA_ARGS) $(MPQCRUN_ARGS) $(CHECK2INPUTS:$(INP)/%=%)

check2_chk:
	$(CHECKOUT) -r $(SRCDIR)/ref $(CHECK2OUTPUTS:$(INP)/%=$(RUN)/%)

.PHONY: inputs
inputs:: h2o h2omp2 mp2r12 h2ofrq basis1 basis2 opt optts symm1 symm2 symm3 ckpt mbpt
inputs:: methods clscf uscf hsosscf input dft orthog ccaint

.PHONY: h2o
h2o: stamp-h2o

.PHONY: h2omp2
h2omp2: stamp-h2omp2

.PHONY: mp2r12
mp2r12: stamp-mp2r12

.PHONY: h2ofrq
h2ofrq: stamp-h2ofrq

.PHONY: orthog
orthog: stamp-orthog

.PHONY: basis1
basis1: stamp-basis1

.PHONY: dft
dft: stamp-dft

.PHONY: basis2
basis2: stamp-basis2

.PHONY: symm1
symm1: stamp-symm1

.PHONY: symm2
symm2: stamp-symm2

.PHONY: symm3
symm3: stamp-symm3

.PHONY: clscf
clscf: stamp-clscf

.PHONY: uscf
uscf: stamp-uscf

.PHONY: hsosscf
hsosscf: stamp-hsosscf

.PHONY: opt
opt: stamp-opt

.PHONY: optts
optts: stamp-optts

.PHONY: ccaint
ccaint: stamp-ccaint

.PHONY: ckpt
ckpt: $(CKPTINPUTS) $(CKPTINPUTS:%.in=%.qci)

.PHONY: mbpt
mbpt: $(MBPTINPUTS) $(MBPTINPUTS:%.in=%.qci)

.PHONY: methods
methods: $(METHODSINPUTS) $(METHODSINPUTS:%.in=%.qci)

.PHONY: input
input: $(INPUTINPUTS) $(INPUTINPUTS:%.in=%.qci)

stamp-h2o: $(INP) $(H2OMASTER) $(INPUTGENDEPS)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(INP) $(H2OMASTER)
	touch stamp-h2o

stamp-h2omp2: $(INP) $(H2OMP2MASTER) $(INPUTGENDEPS)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(INP) $(H2OMP2MASTER)
	touch stamp-h2omp2

stamp-mp2r12: $(INP) $(MP2R12MASTER) $(INPUTGENDEPS)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(INP) $(MP2R12MASTER)
	touch stamp-mp2r12

stamp-h2ofrq: $(INP) $(H2OFRQMASTER) $(INPUTGENDEPS)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(INP) $(H2OFRQMASTER)
	touch stamp-h2ofrq

stamp-orthog: $(INP) $(ORTHOGMASTER) $(INPUTGENDEPS)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(INP) $(ORTHOGMASTER)
	touch stamp-orthog

stamp-basis1: $(INP) $(BASIS1MASTER) $(INPUTGENDEPS)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(INP) $(BASIS1MASTER)
	touch stamp-basis1

stamp-dft: $(INP) $(DFTMASTER) $(INPUTGENDEPS)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(INP) $(DFTMASTER)
	touch stamp-dft

stamp-basis2: $(INP) $(BASIS2MASTER) $(INPUTGENDEPS)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(INP) $(BASIS2MASTER)
	touch stamp-basis2

stamp-opt: $(INP) $(OPTMASTER) $(INPUTGENDEPS)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(INP) $(OPTMASTER)
	touch stamp-opt

stamp-optts: $(INP) $(OPTTSMASTER) $(INPUTGENDEPS)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(INP) $(OPTTSMASTER)
	touch stamp-optts

stamp-symm1: $(INP) $(SYMM1MASTER) $(INPUTGENDEPS)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(INP) $(SYMM1MASTER)
	touch stamp-symm1

stamp-symm2: $(INP) $(SYMM2MASTER) $(INPUTGENDEPS)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(INP) $(SYMM2MASTER)
	touch stamp-symm2

stamp-symm3: $(INP) $(SYMM3MASTER) $(INPUTGENDEPS)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(INP) $(SYMM3MASTER)
	touch stamp-symm3

stamp-clscf: $(INP) $(CLSCFMASTER) $(INPUTGENDEPS)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(INP) $(CLSCFMASTER)
	touch stamp-clscf

stamp-uscf: $(INP) $(USCFMASTER) $(INPUTGENDEPS)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(INP) $(USCFMASTER)
	touch stamp-uscf

stamp-hsosscf: $(INP) $(HSOSSCFMASTER) $(INPUTGENDEPS)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(INP) $(HSOSSCFMASTER)
	touch stamp-hsosscf

stamp-ccaint: $(INP) $(CCAINTMASTER) $(INPUTGENDEPS)
	$(MAKEIN) -p $(PROGRAM) -I$(SRCDIR) -d $(INP) $(CCAINTMASTER)
	touch stamp-ccaint

##############################################################################

.PHONY: checkrun
checkrun: checkckpt checkmbpt checksymm1 checksymm2 checksymm3 checkdft
	@echo Check for complete outputs in the run directory:
	for file in $(ALLOUTPUTS); \
	  do \
	    ($(CHECKOUT) $${file} ); \
	  done

.PHONY: checkh2ofrq
checkh2ofrq:
	@echo Consistency checks on h2o frequency tests:
	@$(CHECKOUT) $(RUN)/h2ofrq_scfsto3g{c1,c2v}frq.out
	@$(CHECKOUT) $(RUN)/h2ofrq_scfsto3g{c1,c2v}optfrq.out
	@$(CHECKOUT) $(RUN)/h2ofrq_mp200sto3g{c1,c2v}frq.out
	@$(CHECKOUT) $(RUN)/h2ofrq_mp200sto3g{c1,c2v}optfrq.out
	@$(CHECKOUT) $(RUN)/h2ofrq_scf6311gss{c1,c2v}frq.out
	@$(CHECKOUT) $(RUN)/h2ofrq_scf6311gss{c1,c2v}optfrq.out
	@$(CHECKOUT) $(RUN)/h2ofrq_mp2006311gss{c1,c2v}frq.out
	@$(CHECKOUT) $(RUN)/h2ofrq_mp2006311gss{c1,c2v}optfrq.out

.PHONY: checkckpt
checkckpt:
	@echo Consistency checks on checkpoint tests:
	@$(CHECKOUT) $(RUN)/ckpt_1hsosscf.out $(RUN)/ckpt_hsosscf.out
	@$(CHECKOUT) $(RUN)/ckpt_1clscf.out $(RUN)/ckpt_clscf.out
	@$(CHECKOUT) $(RUN)/ckpt_1mp2.out $(RUN)/ckpt_mp2.out
	@$(CHECKOUT) $(RUN)/ckpt_qnewtopt.out $(RUN)/ckpt_2qnewtopt.out
	@$(CHECKOUT) $(RUN)/ckpt_efcopt.out $(RUN)/ckpt_2efcopt.out
	@$(CHECKOUT) $(RUN)/ckpt_0zapt2.out $(RUN)/ckpt_1zapt2.out
	@$(CHECKOUT) $(RUN)/ckpt_clksxalpha.out $(RUN)/ckpt_1clksxalpha.out
	@$(CHECKOUT) $(RUN)/ckpt_clksb3lyp.out $(RUN)/ckpt_1clksb3lyp.out
	@$(CHECKOUT) $(RUN)/ckpt_clksbp86.out $(RUN)/ckpt_1clksbp86.out
	@$(CHECKOUT) $(RUN)/ckpt_clkshfg96.out $(RUN)/ckpt_1clkshfg96.out
	@$(CHECKOUT) $(RUN)/ckpt_clksmpwpw91.out $(RUN)/ckpt_1clksmpwpw91.out
	@$(CHECKOUT) $(RUN)/ckpt_clkspbe.out $(RUN)/ckpt_1clkspbe.out
	@$(CHECKOUT) $(RUN)/ckpt_clkspw91.out $(RUN)/ckpt_1clkspw91.out
	@$(CHECKOUT) $(RUN)/ckpt_clksspz81.out $(RUN)/ckpt_1clksspz81.out
	@$(CHECKOUT) $(RUN)/ckpt_clkssvwn1.out $(RUN)/ckpt_1clkssvwn1.out
	@$(CHECKOUT) $(RUN)/ckpt_clkssvwn1rpa.out $(RUN)/ckpt_1clkssvwn1rpa.out
	@$(CHECKOUT) $(RUN)/ckpt_clkssvwn2.out $(RUN)/ckpt_1clkssvwn2.out
	@$(CHECKOUT) $(RUN)/ckpt_clkssvwn3.out $(RUN)/ckpt_1clkssvwn3.out
	@$(CHECKOUT) $(RUN)/ckpt_clkssvwn4.out $(RUN)/ckpt_1clkssvwn4.out
	@$(CHECKOUT) $(RUN)/ckpt_clkssvwn5.out $(RUN)/ckpt_1clkssvwn5.out
	@$(CHECKOUT) $(RUN)/ckpt_hsosksxalpha.out $(RUN)/ckpt_1hsosksxalpha.out
	@$(CHECKOUT) $(RUN)/ckpt_uksxalpha.out $(RUN)/ckpt_1uksxalpha.out

.PHONY: checkmbpt
checkmbpt:
	@echo Consistency checks on MBPT tests:
	@$(CHECKOUT) $(RUN)/mbpt_zapt2{v1,v2}.out
	@$(CHECKOUT) $(RUN)/mbpt_zapt2{v1,v2lb}.out
	@$(CHECKOUT) $(RUN)/mbpt_opt12{v1,v2}.out
	@$(CHECKOUT) $(RUN)/mbpt_opt12{v1,v2lb}.out
	@$(CHECKOUT) $(RUN)/mbpt_opt22{v1,v2}.out
	@$(CHECKOUT) $(RUN)/mbpt_opt22{v1,v2lb}.out
	@$(CHECKOUT) $(RUN)/mbpt_mp2{v1,v2}.out
	@$(CHECKOUT) $(RUN)/mbpt_mp2{v1,v2lb}.out
	@$(CHECKOUT) $(RUN)/mbpt_mp2{v1,mem}.out
	@$(CHECKOUT) $(RUN)/mbpt_mp2{v1,mem_c1}.out
	@$(CHECKOUT) $(RUN)/mbpt_mp2{v1,mem_auto}.out
	@$(CHECKOUT) $(RUN)/mbpt_mp2r12_ne2{,_posix}.out
	@$(CHECKOUT) $(RUN)/mbpt_mp2r12_ne2{,_multipass}.out

.PHONY: checkdft
checkdft:
	@echo Consistency checks on DFT tests:
	for mol in h2 lih beh2 b2h6 nh3 ch4 c2h4 c2h2 h2o hf \
                   nah mgh2 alh sih2 ph3 h2s hcl; \
	do (\
	  $(CHECKOUT) $(RUN)/dft_$${mol}hfs{ultrafine,}631gs*.out \
	); done

.PHONY: checksymm1
checksymm1:
	@echo Consistency checks on symmetry test series 1:
	@for b in sto3g; \
	do (\
	  $(CHECKOUT) $(RUN)/symm1_{c2h2scf$${b}c1,c2h2scf$${b}ci}.out;\
	  $(CHECKOUT) $(RUN)/symm1_{c2h2scf$${b}c1,c2h2scf$${b}c2}.out;\
	  $(CHECKOUT) $(RUN)/symm1_{c2h2scf$${b}c1,c2h2scf$${b}cs}.out;\
	  $(CHECKOUT) $(RUN)/symm1_{c2h2scf$${b}c1,c2h2scf$${b}c2v}.out;\
	  $(CHECKOUT) $(RUN)/symm1_{c2h2scf$${b}c1,c2h2scf$${b}d2h}.out;\
	  $(CHECKOUT) $(RUN)/symm1_{cubscf$${b}c1,cubscf$${b}ci}.out;\
	  $(CHECKOUT) $(RUN)/symm1_{cubscf$${b}c1,cubscf$${b}c2}.out;\
	  $(CHECKOUT) $(RUN)/symm1_{cubscf$${b}c1,cubscf$${b}cs}.out;\
	  $(CHECKOUT) $(RUN)/symm1_{cubscf$${b}c1,cubscf$${b}c2v}.out;\
	  $(CHECKOUT) $(RUN)/symm1_{cubscf$${b}c1,cubscf$${b}d2h}.out;\
	  $(CHECKOUT) $(RUN)/symm1_{c2h2mp222$${b}c1,c2h2mp222$${b}ci}.out;\
	  $(CHECKOUT) $(RUN)/symm1_{c2h2mp222$${b}c1,c2h2mp222$${b}c2}.out;\
	  $(CHECKOUT) $(RUN)/symm1_{c2h2mp222$${b}c1,c2h2mp222$${b}cs}.out;\
	  $(CHECKOUT) $(RUN)/symm1_{c2h2mp222$${b}c1,c2h2mp222$${b}c2v}.out;\
	  $(CHECKOUT) $(RUN)/symm1_{c2h2mp222$${b}c1,c2h2mp222$${b}d2h}.out;\
	  $(CHECKOUT) $(RUN)/symm1_{cubmp284$${b}c1,cubmp284$${b}ci}.out;\
	  $(CHECKOUT) $(RUN)/symm1_{cubmp284$${b}c1,cubmp284$${b}c2}.out;\
	  $(CHECKOUT) $(RUN)/symm1_{cubmp284$${b}c1,cubmp284$${b}cs}.out;\
	  $(CHECKOUT) $(RUN)/symm1_{cubmp284$${b}c1,cubmp284$${b}c2v}.out;\
	  $(CHECKOUT) $(RUN)/symm1_{cubmp284$${b}c1,cubmp284$${b}d2h}.out;\
	); done

.PHONY: checksymm2
checksymm2:
	@echo Consistency checks on symmetry test series 2:
	@for b in ccpvdz ccpvtz ccpvqz ccpv5z; \
	do (\
	  $(CHECKOUT) $(RUN)/symm2_{c2h2scf$${b}c1,c2h2scf$${b}ci}.out;\
	  $(CHECKOUT) $(RUN)/symm2_{c2h2scf$${b}c1,c2h2scf$${b}c2}.out;\
	  $(CHECKOUT) $(RUN)/symm2_{c2h2scf$${b}c1,c2h2scf$${b}cs}.out;\
	  $(CHECKOUT) $(RUN)/symm2_{c2h2scf$${b}c1,c2h2scf$${b}c2v}.out;\
	  $(CHECKOUT) $(RUN)/symm2_{c2h2scf$${b}c1,c2h2scf$${b}d2h}.out;\
	  $(CHECKOUT) $(RUN)/symm2_{c2h2nsscf$${b}auto,c2h2scf$${b}c1}.out;\
	  $(CHECKOUT) $(RUN)/symm2_{cubscf$${b}c1,cubscf$${b}ci}.out;\
	  $(CHECKOUT) $(RUN)/symm2_{cubscf$${b}c1,cubscf$${b}c2}.out;\
	  $(CHECKOUT) $(RUN)/symm2_{cubscf$${b}c1,cubscf$${b}cs}.out;\
	  $(CHECKOUT) $(RUN)/symm2_{cubscf$${b}c1,cubscf$${b}c2v}.out;\
	  $(CHECKOUT) $(RUN)/symm2_{cubscf$${b}c1,cubscf$${b}d2h}.out;\
	); done

.PHONY: checksymm3
checksymm3:
	@for o in $(SYMM3OUTPUTS); \
	do ( \
	  expect=`echo $${o}|sed 's/^.*symm3_.*_\(.*\)_.*$$/\1/'`; \
	  actual=`grep 'g point g' $${o}|sed 's/^.*to \(.*\)$$/\1/'`; \
	  echo $${o}: expected $${expect} got $${actual}; \
	); done

# check all outputs in the ref directory with outputs in the run directory
.PHONY: check
check:
	@echo Comparing outputs in run and ref directories:
	$(CHECKOUT) -d $(REF) $(RUN)

# Check all outputs files defined in this makefile in the ref directory
# with those in the run directory.  This will break on some platforms
# since a very long argument list is used.
.PHONY: fastcheck
checkall:
	@echo Comparing outputs in run and ref directories:
	$(CHECKOUT) -r $(SRCDIR)/ref $(ALLOUTPUTS)

.PHONY: diff
diff:
	-diff -ur $(REF) $(RUN)

##############################################################################

$(RUN):
	mkdir -p $(RUN)

$(INP):
	mkdir -p $(INP)

$(INP)/ckpt_%.in: ckpt/ckpt_%.in
	@mkdir -p $(INP)
	/bin/cp $< $@

$(INP)/ckpt_0clks%.qci: ckpt/energy.qci
	@mkdir -p $(INP)
	/bin/cp $< $@

$(INP)/ckpt_1clks%.qci: ckpt/gradient.qci
	@mkdir -p $(INP)
	/bin/cp $< $@

$(INP)/ckpt_clks%.qci: ckpt/gradient.qci
	@mkdir -p $(INP)
	/bin/cp $< $@

$(INP)/ckpt_0hsosks%.qci: ckpt/energy.qci
	@mkdir -p $(INP)
	/bin/cp $< $@

$(INP)/ckpt_1hsosks%.qci: ckpt/gradient.qci
	@mkdir -p $(INP)
	/bin/cp $< $@

$(INP)/ckpt_hsosks%.qci: ckpt/gradient.qci
	@mkdir -p $(INP)
	/bin/cp $< $@

$(INP)/ckpt_0uks%.qci: ckpt/energy.qci
	@mkdir -p $(INP)
	/bin/cp $< $@

$(INP)/ckpt_1uks%.qci: ckpt/gradient.qci
	@mkdir -p $(INP)
	/bin/cp $< $@

$(INP)/ckpt_uks%.qci: ckpt/gradient.qci
	@mkdir -p $(INP)
	/bin/cp $< $@

$(INP)/ckpt_%.qci: ckpt/ckpt_%.qci
	@mkdir -p $(INP)
	/bin/cp $< $@

$(INP)/mbpt_%.in: mbpt/mbpt_%.in
	@mkdir -p $(INP)
	/bin/cp $< $@

$(INP)/mbpt_%.qci: mbpt/mbpt_%.qci
	@mkdir -p $(INP)
	/bin/cp $< $@

$(INP)/methods_%.in: methods/%.in
	@mkdir -p $(INP)
	/bin/cp $< $@

$(INP)/methods_%.qci: methods/%.qci
	@mkdir -p $(INP)
	/bin/cp $< $@

$(INP)/input_%.in: input/%.in
	@mkdir -p $(INP)
	/bin/cp $< $@

$(INP)/input_%opt.qci:
	@mkdir -p $(INP)
	(echo "method: generic"; echo "optimize: yes") > $@

$(INP)/input_%.qci:
	@mkdir -p $(INP)
	(echo "method: generic"; echo "optimize: no") > $@

##############################################################################

.PHONY: clean
clean: check_clean
	/bin/rm stamp-*
