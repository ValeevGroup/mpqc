
TOPDIR=../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
INCLUDE += -I. -I$(SRCDIR)
CXXINCLUDE += -I. -I$(SRCDIR)

TARGET_TO_MAKE = mkclasses
BIN_OR_LIB = BIN
LD = $(CXX)

DISTFILES = $(CXXFILES)

CXXSRC = main.cc mkclasses.cc db.cc file.cc class.cc \
         ctor.cc members.cc keyval.cc state.cc parent.cc options.cc \
         datamem.cc

GENCXXSRC =parse.cc scan.cc

BINOBJ = $(CXXSRC:%.cc=%.o) $(GENCXXSRC:%.cc=%.o)

INC = mkclasses.h
GENINC = parse.h

DEPENDINCLUDE = $(INC) $(GENINC)

default:: mkclasses

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(BINOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.o=.d)
endif

test: tclasses.h tclasses.cc
	cat tclasses.h
	cat tclasses.cc

tclasses.h tclasses.cc:  tclasses.kv ./mkclasses
	./mkclasses --gensrc --geninc $<

# only works with bison and flex
parse.cc parse.h: parse.yy
	$(YACC) -v -d -o parse.cc $^
	/bin/mv parse.cc.h parse.h
	-@rm -f parse.cc.output

scan.cc: scan.ll
	$(LEX) -L -t $^ | grep -v "extern FILE .yyin" \
	                | grep -v "static int yy_get_next_buffer.*;" \
	                | sed "s/static int yy_get_next_buffer/int yy_get_next_buffer/" \
	                | grep -v "static void yyunput.*;" \
	                | sed "s/static void yyunput/void yyunput/" \
	                | grep -v "static int yyinput.*;" \
	                | sed "s/static int yyinput/int yyinput/" \
	                > scan.cc

