TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

ifeq ($(YACC),bison)
CXXFLAGS += -DBISON
CXXFLAGS_DEBUG += -DBISON
endif
ifeq ($(LEX),flex)
CXXFLAGS += -DFLEX
CXXFLAGS_DEBUG += -DFLEX
endif

MAPIMPL = RAVL

CXXSRC = ipv2.cc ipv2_alloc.cc ipv2_cwk.cc ipv2_data.cc ipv2_error.cc \
	 ipv2_karray.cc ipv2_print.cc ipv2_read.cc \
	 keyval.cc keyvalipv2.cc keyvalval.cc keyvalkey.cc keyvalass.cc \
         keyvalstr.cc keyvalagg.cc keyvalpre.cc

GENCXXSRC =ipv2_parse.cc ipv2_scan.cc \
           keyval$(MAPIMPL)Map.cc keyvalMap.cc

LIBOBJ = $(CXXSRC:%.cc=%.o) $(GENCXXSRC:%.cc=%.o)
GENSRC = $(GENCXXSRC)

INC = ipv2.h keyval.h keyvaldefs.h
GENINC = ipv2_parse.h \
         keyvalImplMap.h keyval$(MAPIMPL)Map.h keyvalMap.h


DEPENDINCLUDE = $(INC) $(GENINC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libkeyval

DISTFILES = Makefile $(CXXSRC) $(INC) ipv2_parse.yy ipv2_scan.ll LIBS.h

default:: $(DEPENDINCLUDE)
default:: libkeyval.a

keyvaltest: keyvaltest.o libkeyval.a ../container/libcontainer.a ../class/libclass.a
	$(CXX) $(CXXFLAGS) -o keyvaltest $^

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
OTHEROBJ = keyvaltest.o
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(OTHEROBJ:.o=.d)
endif

#### yacc and lex ####
# (only works with bison and flex)
ipv2_parse.cc ipv2_parse.h: ipv2_parse.yy
	$(YACC) -v -d -o ipv2_parse.cc $^
	/bin/mv ipv2_parse.cc.h ipv2_parse.h
	-@rm -f ipv2_parse.cc.output

ipv2_scan.cc: ipv2_scan.ll
	$(LEX) -L -t $^ | grep -v "extern FILE .yyin" > ipv2_scan.cc

#############################


keyvalMap.cc keyvalMap.h:
	genclass -2 KeyValKeyword val RefKeyValValue val Map keyval
	prepend-header keyvalMap.h "#include \"keyval.h\""

keyval$(MAPIMPL)Map.cc keyval$(MAPIMPL)Map.h:
	genclass -2 KeyValKeyword val RefKeyValValue val $(MAPIMPL)Map keyval
	prepend-header keyval$(MAPIMPL)Map.h "#include \"keyval.h\""

keyvalImplMap.h:
	echo "#include "\"keyval$(MAPIMPL)Map.h\" > keyvalImplMap.h
	echo "#define MAPCTOR KeyValKeywordRefKeyValValue"$(MAPIMPL)"Map(0)" >> keyvalImplMap.h
