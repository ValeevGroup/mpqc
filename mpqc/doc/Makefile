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
       idxcls.texi idxmac.texi ref.texi state.texi symmetry.texi

# These are the files generated from keyval.h
KEYVALDOC = KeyVal.texi KeyVal.pub.texi KeyVal.pro.texi

# from math/scmat/matrix.h:
MATRIXDOC = RefSCMatrix.texi RefSCMatrix.pub.texi \
            RefSymmSCMatrix.texi RefSymmSCMatrix.pub.texi \
            RefDiagSCMatrix.texi RefDiagSCMatrix.pub.texi \
            RefSCVector.texi RefSCVector.pub.texi \
            RefSCDimension.texi RefSCDimension.pub.texi
#            RefLocalSCMatrix.texi RefLocalSCMatrix.pub.texi \
#            RefSCMatrixKit.texi RefSCMatrixKit.pub.texi

# from math/scmat/abstract.h:
ABSTRACTDOC = SCMatrix.texi SCMatrix.pub.texi \
              SymmSCMatrix.texi SymmSCMatrix.pub.texi \
              DiagSCMatrix.texi DiagSCMatrix.pub.texi \
              SCVector.texi SCVector.pub.texi \
              SCDimension.texi SCDimension.pub.texi \
              SCMatrixKit.texi SCMatrixKit.pub.texi
#              LocalSCMatrix.texi LocalSCMatrix.pub.texi \

# from math/symmetry/pointgroup.h
SYMMETRYDOC = SymmetryOperation.texi SymmetryOperation.pub.texi \
              IrreducibleRepresentation.texi \
              IrreducibleRepresentation.pub.texi \
              IrreducibleRepresentation.pri.texi \
              CharacterTable.texi CharacterTable.pub.texi \
              CharacterTable.pri.texi \
              PointGroup.texi PointGroup.pub.texi

# from chemistry/molecule/chemelem.h
CHEMEDOC = ChemicalElement.texi ChemicalElement.pub.texi

# from chemistry/molecule/chemelem.h
ATOMDOC = AtomicCenter.texi AtomicCenter.pub.texi

# from chemistry/molecule/molecule.h
MOLDOC = Molecule.texi Molecule.pub.texi Molecule.pri.texi

# from chemistry/molecule/coor.h
COORDOC = IntCoor.texi IntCoor.pub.texi \
          SumIntCoor.texi SumIntCoor.pub.texi \
          SetIntCoor.texi SetIntCoor.pub.texi

# These are all of texi files used to generate the documentation.
# The top level document, SC.texi, must be first in the list.
GENDOC = $(KEYVALDOC) $(MATRIXDOC) $(ABSTRACTDOC) $(SYMMETRYDOC) \
         $(CHEMEDOC) $(ATOMDOC) $(MOLDOC) $(COORDOC)

ALLDOC = $(DOC) $(GENDOC)

TEXIDIRS = $(SRCDIR)
EXTRATEXINPUT = $(shell /bin/sh -f $(SRCDIR)/addcolon.sh $(TEXIDIRS))
INCLUDEDIRS = $(shell /bin/sh -f $(SRCDIR)/addi.sh $(TEXIDIRS))

default: info

.PHONY:	info html dvisee pssee SC.dvi

info: $(ALLDOC)
	makeinfo $(INCLUDEDIRS) $<

html: $(ALLDOC)
	#TEXINPUTS=$(EXTRATEXINPUT):$$TEXINPUTS & perl $(T2H) -verbose -deepsplit -menu $<
	TEXINPUTS=$(EXTRATEXINPUT):$$TEXINPUTS & perl $(T2H) -verbose -split_chapter -split_node -menu $<

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
	/bin/rm -f SC.aux SC.cl SC.cls SC.cp SC.cps SC.dvi SC.fn SC.ky SC.log
	/bin/rm -f SC.mc SC.mcs
	/bin/rm -f SC.pg SC.toc SC.tp SC.vr SC.ps
	/bin/rm -f SC.*.bak
	/bin/rm -f $(GENDOC)

$(KEYVALDOC): $(SRCLIBDIR)/util/keyval/keyval.h
	$(C2T) $<

$(MATRIXDOC): $(SRCLIBDIR)/math/scmat/matrix.h
	$(C2T) $<

$(SYMMETRYDOC): $(SRCLIBDIR)/math/symmetry/pointgrp.h
	$(C2T) $<

$(ABSTRACTDOC): $(SRCLIBDIR)/math/scmat/abstract.h
	$(C2T) $<

$(CHEMEDOC): $(SRCLIBDIR)/chemistry/molecule/chemelem.h
	$(C2T) $<

$(ATOMDOC): $(SRCLIBDIR)/chemistry/molecule/atomcent.h
	$(C2T) $<

$(MOLDOC): $(SRCLIBDIR)/chemistry/molecule/molecule.h
	$(C2T) $<

$(COORDOC): $(SRCLIBDIR)/chemistry/molecule/coor.h
	$(C2T) $<

