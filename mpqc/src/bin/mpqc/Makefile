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

#################################################################

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

mpqc: $(BINOBJ) $(LIBS)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(SYSLIBS) $(LTLINKBINOPTS)

$(TOPDIR)/lib/libmpqc.a: $(BINOBJ)
	/bin/rm -f $@
	$(AR) $(ARFLAGS) $@ $^
	$(RANLIB) $@

$(TOPDIR)/lib/libmpqc.la: $(BINOBJ)
	$(LTLINK) $(CXX) -o $@ $^ $(LTLINKLIBOPTS)

install_devel:: $(TOPDIR)/lib/libmpqc.$(LIBSUF)
	$(INSTALL) $(INSTALLDIROPT) $(libdir)
	$(LTINST) $(INSTALL) $(INSTALLLIBOPT) $< $(libdir)

install:: mpqc
	$(INSTALL) $(INSTALLDIROPT) $(bindir)
	$(LTINST) $(INSTALL) $(INSTALLBINOPT) $< $(bindir)

install::
	sed "s+exec perl+exec perl -I$(prefix)/lib/perl+" \
		< $(SRCDIR)/validate/checkout.pl \
		> $(bindir)/chkmpqcout
	chmod 555 $(bindir)/chkmpqcout

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
	echo "#include <iostream>" >> scan.cc
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
	                | grep -v "extern.*isatty" \
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
