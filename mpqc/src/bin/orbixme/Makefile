TOPDIR=../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
INCLUDE += -I. -I$(SRCDIR)
CXXINCLUDE += -I. -I$(SRCDIR)

IDL = idl -h .h -s serv.cc -c clnt.cc

ITCLT           = -lITclt
IRCLT           = -lIRclt
ITSRV           = -lITsrv

ORBDIR = /net/herzberg/usr/people/ejfried/Orbix1.3
ORBINCLUDE = -I$(ORBDIR)/include
ORBLIBDIR = -L$(ORBDIR)/lib

CXXINCLUDE += $(ORBINCLUDE)

ORBSERVERLIBS = $(ORBLIBDIR) $(ITSRV)
ORBCLIENTLIBS = $(ORBLIBDIR) $(ITCLT)

SERVEROBJ = server.o moleserv.o molei.o
CLIENTOBJ = client.o moleclnt.o

INC = molei.h
GENINC = mole.h
GENSRC = moleserv.cc moleclnt.cc

DEPENDINCLUDE = $(INC) $(GENINC)

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

.PHONY: default
default:: moleserv moleclnt

moleserv: $(SERVEROBJ) $(LIBS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(SYSLIBS) $(ORBSERVERLIBS)

moleclnt: $(CLIENTOBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(ORBCLIENTLIBS)

%.o: %.C
	$(CXX) $(CXXFLAGS) -c -o $@ $<

moleclnt.cc moleserv.cc mole.h: mole.idl
	$(IDL) -B $^

oclean::
	/bin/rm -f moleserv moleclnt

distclean::
	/bin/rm -f $(GENINC) $(GENSRC)


#################################################################

default:: $(DEPENDINCLUDE)

$(CLIENTOBJ:.o=.d) $(SERVEROBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(CLIENTOBJ:.o=.d) $(BINOBJ:.o=.d)
endif
