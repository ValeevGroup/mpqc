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
XLIBS =
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_PSI),yes)
XDEF += -DHAVE_CHEMISTRY_QC_PSI
endif
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_CC),yes)
XDEF += -DHAVE_CHEMISTRY_QC_CC
endif
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_CINTS),yes)
 ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_QC_MBPTR12),yes)
 XDEF += -DHAVE_CHEMISTRY_QC_MBPTR12
 endif
XDEF += -DHAVE_CHEMISTRY_QC_CINTS
endif

ifeq ($(HAVE_LIBINT),yes)
  XLIBS += -lint
endif
ifeq ($(HAVE_LIBR12),yes)
  XLIBS += -lr12
endif
ifeq ($(HAVE_LIBDERIV),yes)
  XLIBS += -lderiv
endif

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(XDEF) $(SRCDIR)/LIBS.h)


#################################################################

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

mpqc: $(BINOBJ) $(LIBS)
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(XLIBS) $(SYSLIBS) $(LTLINKBINOPTS)

$(TOPDIR)/lib/libmpqc.a: $(BINOBJ)
	/bin/rm -f $@
	$(AR) $(ARFLAGS) $@ $^
	$(RANLIB) $@

$(TOPDIR)/lib/libmpqc.la: $(BINOBJ)
	$(LTLINK) $(CXX) -o $@ $^ $(LTLINKLIBOPTS)

install_devel:: $(TOPDIR)/lib/libmpqc.$(LIBSUF)
	$(INSTALL) $(INSTALLDIROPT) $(installroot)$(libdir)
	$(LTINST) $(INSTALL) $(INSTALLLIBOPT) $< $(installroot)$(libdir)

install:: mpqc
	$(INSTALL) $(INSTALLDIROPT) $(installroot)$(bindir)
	$(LTINST) $(INSTALL) $(INSTALLBINOPT) $< $(installroot)$(bindir)

install::
	/bin/rm -f $(installroot)$(bindir)/chkmpqcout
	sed "s+exec perl+exec perl -I$(scdatadir)/perl+" \
		< $(SRCDIR)/validate/checkout.pl \
		> $(installroot)$(bindir)/chkmpqcout
	chmod 555 $(installroot)$(bindir)/chkmpqcout

clean::
	/bin/rm -f mpqc
	/bin/rm -f parse.tmp*

distclean::
	/bin/rm -f mpqc
	/bin/rm -f parse.tmp*
	/bin/rm -f $(GENCXXSRC) $(GENINC) parse.output

# (only works with bison and flex)
parse.cc: parse.yy
	$(BISON) -v -d -o parse.tmp.cc $^
	cat parse.tmp.cc \
		| sed "s/^int yyparse.*;$$//" \
		| sed "s/^YYPARSE_RETURN_TYPE yyparse.*;$$//" \
		| sed "s/yylval/MPQCInylval/g" \
		> parse.cc
	if test -f parse.tmp.cc.h; then \
	  echo "Older bison detected."; \
	  cat parse.tmp.cc.h | sed "s/yylval/MPQCInylval/g" > parse.h; \
	  /bin/mv parse.tmp.cc.output parse.output; \
	else \
	  echo "Newer bison detected."; \
	  cat parse.tmp.hh | sed "s/yylval/MPQCInylval/g" > parse.h; \
	  /bin/mv parse.tmp.output parse.output; \
	fi
	#-@rm -f parse.tmp*

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
