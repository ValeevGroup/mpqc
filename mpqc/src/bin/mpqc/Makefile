TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
include $(TOPDIR)/lib/Makedirlist
INCLUDE += -I. -I$(SRCDIR)
CXXINCLUDE += -I. -I$(SRCDIR)

TARGET_TO_MAKE = mpqc

CXXSRC = mpqc.cc mpqcin.cc
GENCXXSRC = parse.cc scan.cc
CSRC =
INC =
GENINC = parse.h

DEPENDINCLUDE = $(INC)

BINOBJ = $(CXXSRC:%.cc=%.$(OBJSUF)) $(GENCXXSRC:%.cc=%.$(OBJSUF)) $(CSRC:%.c=%.$(OBJSUF))

DISTFILES = $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h

default:: mpqc

XDEF =
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_PSI),yes)
XDEF += -DHAVE_CHEMISTRY_QC_PSI
endif
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_CC),yes)
XDEF += -DHAVE_CHEMISTRY_QC_CC
endif

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(XDEF) $(SRCDIR)/LIBS.h)

ifeq ($(ARCH),I860)
  DEFINES += -DLINPACKD=-lkmath
endif

#################################################################

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

mpqc: $(BINOBJ) $(LIBS)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

install:: mpqc
	$(INSTALL) $(INSTALLDIROPT) $(bindir)
	$(LTINST) $(INSTALL) $(INSTALLBINOPT) $< $(bindir)

distclean::
	/bin/rm -f mpqc
	/bin/rm -f $(GENCXXSRC) $(GENINC) parse.output

# (only works with bison and flex)
parse.cc: parse.yy
	$(BISON) -v -d -o parse.tmp $^
	cat parse.tmp \
		| sed "s/^int yyparse.*;$$//" \
		| sed "s/^YYPARSE_RETURN_TYPE yyparse.*;$$//" \
		| sed "s/yylval/MPQCInylval/g" \
		> parse.cc
	cat parse.tmp.h \
		| sed "s/yylval/MPQCInylval/g" \
		> parse.h
	/bin/mv parse.tmp.output parse.output
	-@rm -f parse.tmp parse.tmp.h

scan.cc: scan.ll
	echo "#ifdef HAVE_CONFIG_H"  > scan.cc
	echo "#include <scconfig.h>" >> scan.cc
	echo "#endif" >> scan.cc
	echo "#include <iostream.h>" >> scan.cc
	echo "#ifdef USING_NAMESPACE_STD" >> scan.cc
	echo "using namespace std;" >> scan.cc
	echo "#endif" >> scan.cc
	$(FLEX) -L -t $^ | grep -v "extern FILE .yyin" \
	                | grep -v "static int yy_get_next_buffer.*;" \
	                | grep -v "class istream;" \
	                | sed "s/static int yy_get_next_buffer/int yy_get_next_buffer/" \
	                | grep -v "static void yyunput.*;" \
	                | sed "s/static void yyunput/void yyunput/" \
	                | grep -v "static int yyinput.*;" \
	                | sed "s/static int yyinput/int yyinput/" \
			| sed "s/yylval/MPQCInylval/g" \
	                >> scan.cc

parse.h: parse.cc

scan.d: parse.h

#################################################################

$(BINOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.$(OBJSUF)=.d)
endif
