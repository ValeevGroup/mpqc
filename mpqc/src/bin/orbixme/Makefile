TOPDIR=../..
ifndef SRCDIR
  SRCDIR=.
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

CHEMSERVSRC = chemserv.cc $(ENERGYSERVSRC) $(MOLECULESERVSRC) $(KEYVALSERVSRC)
CHEMCLNTSRC = chemclnt.cc $(ENERGYCLNTSRC) $(MOLECULECLNTSRC) $(KEYVALCLNTSRC)
CHEMSERVOBJ = $(CHEMSERVSRC:%.cc=%.o)
CHEMCLNTOBJ = $(CHEMCLNTSRC:%.cc=%.o)

GENSRC = $(KEYVALSERVGENSRC) $(KEYVALCLNTGENSRC) \
         $(ENERGYSERVGENSRC) $(ENERGYCLNTGENSRC) \
         $(MOLECULESERVGENSRC) $(MOLECULECLNTGENSRC)
GENINC = $(KEYVALGENINC) $(ENERGYGENINC) $(MOLECULEGENINC)
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
