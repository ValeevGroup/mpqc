TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
INCLUDE += -I. -I$(SRCDIR)

BIN_OR_LIB = BIN
TARGET_TO_MAKE = tmpl

DISTFILES = main.c tmpllex.l tmpl.h

DEPENDCSRC = main.c tmpllex.c

BINOBJ = main.o tmpllex.o

default:: tmpl

tmpllex.c: tmpllex.l
	$(FLEX) -t $^ > tmpllex.c

clean::
	/bin/rm -f tmpllex.c

tmplparse.c: tmplparse.y
	yacc -d -v $^
	mv -f y.tab.c tmplparse.c
	mv -f y.tab.h tmplparse.h

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(BINOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.o=.d)
endif

