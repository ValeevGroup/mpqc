TOPDIR=../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

BIN_OR_LIB = BIN
TARGET_TO_MAKE = tmpl

#SYSLIBS := $(SYSLIBS) -ll

DISTFILES = main.c tmpllex.l tmpl.h

DEPENDCSRC = main.c tmpllex.c

BINOBJ = main.o tmpllex.o

default:: tmpl

# tmpl.h is put into the source include directory.  This is
# necessary because tmpl is not compiled for all architectures.
install_inc::
	$(INCINSTALL) tmpl.h $(SRCDIR)/$(TOPDIR)/include

tmpllex.c: tmpllex.l
	lex -t $^ > tmpllex.c

tmplparse.c: tmplparse.y
	yacc -d -v $^
	mv -f y.tab.c tmplparse.c
	mv -f y.tab.h tmplparse.h

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(BINOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.o=.d)
endif

