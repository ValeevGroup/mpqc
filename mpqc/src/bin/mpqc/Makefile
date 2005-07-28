TOPDIR=../../..
ifndef SRCDIR
  BUILDING_IN_SRCDIR=yes
  SRCDIR=$(shell pwd)
endif

LOCALMAKEFILE_OPTIONAL = yes
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

ifeq ($(LOCALMAKEFILE_FOUND),yes)

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
ifeq ($(HAVE_SC_SRC_LIB_CHEMISTRY_CCA),yes)
  ifeq ($(ENABLESHARED),yes)
    LTLINKBINOPTS += -R$(CCAFE_LIB) -R$(CCA_CHEM_LIB)
  else
    XLIBS += -L$(CCA_CHEM_LIB) -L$(CCAFE_LIB) -L$(CCA_SPEC_BABEL_LIB) -L$(BABEL_LIB) \
             -lccachem_cxx_server -lccachem_cxx_client -lccafeCore -lcca -lsidl -lxml2 -lz
  endif
  DEFINES += -DHAVE_CHEMISTRY_CCA -DCCA_PATH=\"$(libdir)/cca\"
  CPPFLAGS += -I../../lib/chemistry/cca
endif

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(XDEF) $(SRCDIR)/LIBS.h)

#################################################################

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

mpqc: $(BINOBJ) $(LIBS)
	echo "$(LISTLIBS) $(INCLUDE) $(DEFINES) $(XDEF) $(SRCDIR)/LIBS.h"
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
	$(INSTALL) $(INSTALLSCRIPTOPT) mpqcrun $(installroot)$(bindir)/mpqcrun
	$(INSTALL) $(INSTALLSCRIPTOPT) ccarun $(installroot)$(bindir)/ccarun
	$(INSTALL) $(INSTALLDIROPT) $(installroot)$(scdatadir)
	$(INSTALL) $(INSTALLSCRIPTOPT) $(SRCDIR)/mpqcrunproc $(installroot)$(scdatadir)/mpqcrunproc
	$(INSTALL) $(INSTALLSCRIPTOPT) $(SRCDIR)/ccarunproc $(installroot)$(scdatadir)/ccarunproc

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
	/bin/rm -f mpqcrun

ifneq ($(BUILDING_IN_SRCDIR),yes)
ifeq ($(wildcard parse.cc),parse.cc)
$(error "parse.cc exists in an object directory.  This is now longer necessary or allowed.  Delete the file to continue.")
endif
ifeq ($(wildcard scan.cc),scan.cc)
$(error "scan.cc exists in an object directory.  This is now longer necessary or allowed.  Delete the file to continue.")
endif
ifeq ($(wildcard parse.h),parse.h)
$(error "parse.h exists in an object directory.  This is now longer necessary or allowed.  Delete the file to continue.")
endif
endif

$(SRCDIR)/parse.cc: $(SRCDIR)/parse.yy
	@echo WARNING: The file $@ is out of date.
	@echo It can be built by running \"make parser DODEPEND=no\" in the source directory.
	@echo You may also get the message on files checked out of CVS, in which case you can touch $@ to stop getting this message.

$(SRCDIR)/parse.h: $(SRCDIR)/parse.yy
	@echo WARNING: The file $@ is out of date.
	@echo It can be built by running \"make parser DODEPEND=no\" in the source directory.
	@echo You may also get the message on files checked out of CVS, in which case you can touch $@ to stop getting this message.

$(SRCDIR)/scan.cc: $(SRCDIR)/scan.ll
	@echo WARNING: The file $@ is out of date.
	@echo It can be built by running \"make scanner DODEPEND=no\" in the source directory.
	@echo You may also get the message on files checked out of CVS, in which case you can touch $@ to stop getting this message.

#################################################################

$(BINOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.$(OBJSUF)=.d)
endif

endif

ifndef FLEX
FLEX=flex
endif

ifndef BISON
BISON=bison
endif

notobjdir_default:
	@echo "Building in an unconfigured source directory."
	@echo "The following make targets are available:"
	@echo "  make FLEX=<flex> scanner"
	@echo "  make BISON=<bison> parser"
	@echo "Be sure to replace FlexLexer.h in the include"
	@echo "directory with the correct version."

# (only works with bison and flex)
.PHONY: parser
parser:
	$(BISON) -v -d -o parse.tmp.cc parse.yy
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

.PHONY: scanner
scanner:
	echo "#ifdef HAVE_CONFIG_H"  > scan.cc
	echo "#include <scconfig.h>" >> scan.cc
	echo "#endif" >> scan.cc
	echo "#include <iostream>" >> scan.cc
	echo "#ifdef USING_NAMESPACE_STD" >> scan.cc
	echo "using namespace std;" >> scan.cc
	echo "#endif" >> scan.cc
	$(FLEX) -L -t scan.ll | grep -v "extern FILE .yyin" \
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
