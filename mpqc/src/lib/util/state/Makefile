TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
LD = $(CXX)

SETIMPL = AVL
MAPIMPL = RAVL

CXXSRC = state.cc state_bin.cc \
         state_ptr.cc state_text.cc state_xdr.cc state_file.cc qc_xdr_i.cc

GENCXXSRC = stateptrSet.cc \
            stateptr$(SETIMPL)Set.cc statenumSet.cc statenum$(SETIMPL)Set.cc \
            classdMap.cc classd$(MAPIMPL)Map.cc

LIBOBJ = $(CXXSRC:%.cc=%.o) $(GENCXXSRC:%.cc=%.o)
GENSRC = $(GENCXXSRC)

INC =  LIBS.h state.h state_debug.h state_net.h state_ptr.h \
       statenumdefs.h stateptrdefs.h qc_xdr.h classddefs.h \
       qc_xdr_i.h stattmpl.h
GENINC = stateptrImplSet.h statenumImplSet.h \
         stateptrSet.h stateptr$(SETIMPL)Set.h statenumSet.h \
         statenum$(SETIMPL)Set.h \
         classdImplMap.h classdMap.h classd$(MAPIMPL)Map.h \
         statmacr.h

DEPENDINCLUDE = $(INC) $(GENINC)

DISTFILES = $(CXXSRC) $(INC) Makefile LIBS.h

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libstate

default:: $(DEPENDINCLUDE)

interface:: $(DEPENDINCLUDE)

statmacr.h: stattmpl.h $(TMPL2MACDEP)
	$(TMPL2MAC) $< > statmacr.h

statetest.o: statetest.cc
	$(CXX) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

statetest: statetest.o libstate.$(LIBSUF) libkeyval.$(LIBSUF) libclass.$(LIBSUF) libcontainer.$(LIBSUF) libref.$(LIBSUF)
	$(LD) $(LDFLAGS) -o statetest $^ $(SYSLIBS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

statetest.d $(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) statetest.d
endif

stateptrSet.cc stateptrSet.h:
	$(GENCLASS) StateDataPtr val Set stateptr
	prepend-header stateptrSet.h "#include <util/state/state_ptr.h>"
	prepend-header stateptrSet.h "#include <util/state/stateptrdefs.h>"
	/bin/mv stateptrSet.h tmp.stateptrSet.h
	grep -v "#include \"stateptrdefs.h\"" < tmp.stateptrSet.h > stateptrSet.h
	/bin/rm tmp.stateptrSet.h

stateptr$(SETIMPL)Set.cc stateptr$(SETIMPL)Set.h:
	$(GENCLASS) StateDataPtr val $(SETIMPL)Set stateptr
	prepend-header stateptr$(SETIMPL)Set.h "#include <util/state/state_ptr.h>"

stateptrImplSet.h:
	echo "#include <util/state/stateptr$(SETIMPL)Set.h>" > stateptrImplSet.h
	echo "#define StateDataPtr_CTOR StateDataPtr"$(SETIMPL)"Set()" >> stateptrImplSet.h

statenumSet.cc statenumSet.h:
	$(GENCLASS) StateDataNum val Set statenum
	prepend-header statenumSet.h "#include <util/state/state_ptr.h>"
	prepend-header statenumSet.h "#include <util/state/statenumdefs.h>"
	/bin/mv statenumSet.h tmp.statenumSet.h
	grep -v "#include \"statenumdefs.h\"" < tmp.statenumSet.h > statenumSet.h
	/bin/rm tmp.statenumSet.h

statenum$(SETIMPL)Set.cc statenum$(SETIMPL)Set.h:
	$(GENCLASS) StateDataNum val $(SETIMPL)Set statenum
	prepend-header statenum$(SETIMPL)Set.h "#include <util/state/state_ptr.h>"

statenumImplSet.h:
	echo "#include <util/state/statenum$(SETIMPL)Set.h>" > statenumImplSet.h
	echo "#define StateDataNum_CTOR StateDataNum"$(SETIMPL)"Set()" >> statenumImplSet.h

##########################################################################

classdMap.cc classdMap.h:
	$(GENCLASS) -2 ClassDescP val int val Map classd
	prepend-header classdMap.h "#include <util/class/class.h>"
	prepend-header classdMap.h "#include <util/state/classddefs.h>"
	/bin/mv classdMap.h tmp.classdMap.h
	grep -v "#include \"classddefs.h\"" < tmp.classdMap.h > classdMap.h
	/bin/rm tmp.classdMap.h

classd$(MAPIMPL)Map.cc classd$(MAPIMPL)Map.h:
	$(GENCLASS) -2 ClassDescP val int val $(MAPIMPL)Map classd
	prepend-header classd$(MAPIMPL)Map.h "#include <util/class/class.h>"

classdImplMap.h:
	echo "#include <util/state/classd$(MAPIMPL)Map.h>" > classdImplMap.h
	echo "#define ClassDescPintMap_CTOR ClassDescPint"$(MAPIMPL)"Map(0)" >> classdImplMap.h

