TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
INCLUDE += -I. -I$(SRCDIR)
CXXINCLUDE += -I. -I$(SRCDIR)

IDL = idl -h .h -s S.cc -c C.cc

ITCLT           = -lITclt
IRCLT           = -lIRclt
ITSRV           = -lITsrv

ORBDIR = /net/herzberg/usr/people/ejfried/Orbix1.3
ORBINCLUDE = -I$(ORBDIR)/include
ORBLIBDIR = -L$(ORBDIR)/lib

CXXINCLUDE += $(ORBINCLUDE)

ORBSERVERLIBS = $(ORBLIBDIR) $(ITSRV)
ORBCLIENTLIBS = $(ORBLIBDIR) $(ITCLT)

INC = energyi.h keyvali.h

KEYVALGENINC = keyval.h
KEYVALSERVGENSRC = keyvalS.cc
KEYVALCLNTGENSRC = keyvalC.cc
KEYVALSERVSRC = keyvali.cc $(KEYVALSERVGENSRC)
KEYVALCLNTSRC = $(KEYVALCLNTGENSRC)

ENERGYGENINC = energy.h
ENERGYSERVGENSRC = energyS.cc
ENERGYCLNTGENSRC = energyC.cc
ENERGYSERVSRC = energyi.cc $(ENERGYSERVGENSRC)
ENERGYCLNTSRC = $(ENERGYCLNTGENSRC)

MOLECULEGENINC = molecule.h
MOLECULESERVGENSRC = moleculeS.cc
MOLECULECLNTGENSRC = moleculeC.cc
MOLECULESERVSRC = moleculei.cc $(MOLECULESERVGENSRC)
MOLECULECLNTSRC = $(MOLECULECLNTGENSRC)

FUNCTIONGENINC = function.h
FUNCTIONSERVGENSRC = functionS.cc
FUNCTIONCLNTGENSRC = functionC.cc
FUNCTIONSERVSRC = functioni.cc $(FUNCTIONSERVGENSRC)
FUNCTIONCLNTSRC = $(FUNCTIONCLNTGENSRC)

OPTIMIZEGENINC = optimize.h
OPTIMIZESERVGENSRC = optimizeS.cc
OPTIMIZECLNTGENSRC = optimizeC.cc
OPTIMIZESERVSRC = optimizei.cc $(OPTIMIZESERVGENSRC)
OPTIMIZECLNTSRC = $(OPTIMIZECLNTGENSRC)

CHEMSERVSRC = chemserv.cc $(ENERGYSERVSRC) $(MOLECULESERVSRC) \
                          $(FUNCTIONSERVSRC) $(OPTIMIZESERVSRC) \
                          $(KEYVALSERVSRC)
CHEMCLNTSRC = chemclnt.cc $(ENERGYCLNTSRC) $(MOLECULECLNTSRC) \
                          $(FUNCTIONCLNTSRC) $(OPTIMIZECLNTSRC) \
                          $(KEYVALCLNTSRC)
CHEMSERVOBJ = $(CHEMSERVSRC:%.cc=%.o)
CHEMCLNTOBJ = $(CHEMCLNTSRC:%.cc=%.o)

GENSRC = $(KEYVALSERVGENSRC) $(KEYVALCLNTGENSRC) \
         $(ENERGYSERVGENSRC) $(ENERGYCLNTGENSRC) \
         $(MOLECULESERVGENSRC) $(MOLECULECLNTGENSRC) \
         $(FUNCTIONSERVGENSRC) $(FUNCTIONCLNTGENSRC) \
         $(OPTIMIZESERVGENSRC) $(OPTIMIZECLNTGENSRC)
GENINC = $(KEYVALGENINC) $(ENERGYGENINC) $(MOLECULEGENINC) \
         $(OPTIMIZEGENINC) $(FUNCTIONGENINC)
DEPENDINCLUDE = $(INC) $(GENINC)

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

.PHONY: default
default:: chemserv chemclnt

chemserv: $(CHEMSERVOBJ) $(LIBS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SYSLIBS) $(ORBSERVERLIBS)

chemclnt: $(CHEMCLNTOBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(ORBCLIENTLIBS)

%.o: %.C
	$(CXX) $(CXXFLAGS) -c -o $@ $<

energyC.cc energyS.cc energy.h: energy.idl
	/bin/rm -f energyC.cc energyS.cc energy.h
	/bin/rm -f energyC.c++ energyS.c++
	$(IDL) $(CXXINCLUDE) -B $^

keyvalC.cc keyvalS.cc keyval.h: keyval.idl
	/bin/rm -f keyvalS.cc keyvalC.cc keyval.h
	/bin/rm -f keyvalS.c++ keyvalC.c++
	$(IDL) $(CXXINCLUDE) -B $^

optimizeC.cc optimizeS.cc optimize.h: optimize.idl
	/bin/rm -f optimizeS.cc optimizeC.cc optimize.h
	/bin/rm -f optimizeS.c++ optimizeC.c++
	$(IDL) $(CXXINCLUDE) -B $^

functionC.cc functionS.cc function.h: function.idl
	/bin/rm -f functionS.cc functionC.cc function.h
	/bin/rm -f functionS.c++ functionC.c++
	$(IDL) $(CXXINCLUDE) -B $^

moleculeC.cc moleculeS.cc molecule.h: molecule.idl
	/bin/rm -f moleculeS.cc moleculeC.cc molecule.h
	/bin/rm -f moleculeS.c++ moleculeC.c++
	$(IDL) $(CXXINCLUDE) -B $^

oclean::
	/bin/rm -f chemserv chemclnt

distclean::
	/bin/rm -f $(GENINC) $(GENSRC)


#################################################################

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

#################################################################

default:: $(DEPENDINCLUDE)

DEPENDOBJ = $(CHEMSERVOBJ) $(CHEMCLNTOBJ)
$(DEPENDOBJ): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(DEPENDOBJ:.o=.d)
endif
