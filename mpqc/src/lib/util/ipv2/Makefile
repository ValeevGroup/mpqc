TOPDIR = ../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libipv2
BIN_OR_LIB = LIB

ifeq ($(patsubst -S%,,$(LEX)),flex)
DEFINES += -DFLEX
endif

ifeq ($(YACC),bison)
DEFINES += -DBISON
endif

CSRC = ip_alloc.c ip_cwk.c ip_data.c ip_error.c ip_karray.c ip_print.c \
       ip_read.c

OTHERSRC = scan.l parse.y alloca.c

GENSRC = scan.c parse.c

ifeq ($(ARCH),NCUBE_V2)
CSRC += alloca.c
endif
ifeq ($(ARCH),NCUBE_V3)
CSRC += alloca.c
endif

INC = ip_error.h ip_global.h ip_lib.h ip_types.h

GENINC = parse.h ip_alloc.gbl ip_cwk.gbl ip_data.gbl ip_error.gbl \
         ip_karray.gbl ip_print.gbl ip_read.gbl scan.gbl

LIBOBJ = $(CSRC:%.c=%.$(OBJSUF)) $(GENSRC:%.c=%.$(OBJSUF))

DISTFILES = $(CSRC) $(OTHERSRC) $(INC) Makefile LIBS.h

DEPENDINCLUDE = $(INC) $(GENINC)

###################################################################

default:: $(DEPENDINCLUDE)
default:: libipv2.a

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE) $(GENSRC)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif


###################################################################

ifeq ($(YACC),yacc)
parse.c: parse.y
	$(YACC) -d $<
	/bin/mv y.tab.c parse.c
	/bin/mv y.tab.h parse.h

parse.h: parse.c
else
parse.c: parse.y
	$(YACC) $< -d -o parse.c

parse.h: parse.c
endif

scan.$(OBJSUF): scan.c
	$(CC) $(CFLAGS) -c scan.c

scan.c: scan.l parse.h
	$(LEX) -t $< > scan.c
