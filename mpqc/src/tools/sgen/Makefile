TOPDIR=../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = sgen
BIN_OR_LIB = BIN

ifeq ($(LEX),flex)
DEFINES := $(DEFINES) -DFLEX
endif

COPTIONS = $(COPTIONS_DEBUG)
CFLAGS = $(COPTIONS) $(INCLUDE) $(DEFINES)

ifeq ($(ARCH),SGI)
#CC = $(GCC)
#COPTIONS = $(GCOPTIONS)
endif

ifeq ($(LEX),lex)
LIBS = -ly
endif

CFILES = \
alloca.c \
xalloc.c \
send0.c \
recv0.c \
assign.c \
iseq.c \
zero.c \
fread.c \
fwrite.c \
alloc.c \
bcast0.c \
rbcast0.c \
sbcast0.c \
bread.c \
bwrite.c \
gen_read.c \
gen_write.c \
error.c \
free.c \
init.c \
ip.c \
print.c \
sgen.c \
sgen_read.c \
sgen_util.c \
struct.c

LFILES = scan.l
YFILES = parse.y

DISTFILES = $(CFILES) $(LFILES) $(YFILES)

SRC=$(CFILES) $(YFILES:%.y=%.c) $(LFILES:%.l=%.c)

DEPENDCSRC = $(SRC)

BINOBJ = $(SRC:%.c=%.o)

DEPENDINCLUDE =  parse.h \
alloc.gbl alloc.lcl \
assign.gbl assign.lcl \
bcast0.gbl bcast0.lcl \
bread.gbl bread.lcl \
bwrite.gbl bwrite.lcl \
error.gbl error.lcl \
fread.gbl fread.lcl \
free.gbl free.lcl \
fwrite.gbl fwrite.lcl \
gen_read.gbl gen_read.lcl \
gen_write.gbl gen_write.lcl \
init.gbl init.lcl \
ip.gbl ip.lcl \
iseq.gbl iseq.lcl \
print.gbl print.lcl \
rbcast0.gbl rbcast0.lcl \
recv0.gbl recv0.lcl \
sbcast0.gbl sbcast0.lcl \
send0.gbl send0.lcl \
sgen.gbl sgen.lcl \
sgen_read.gbl sgen_read.lcl \
sgen_util.gbl sgen_util.lcl \
sizeof.gbl sizeof.lcl \
struct.gbl struct.lcl \
zero.gbl zero.lcl


default:: sgen

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

scan.o: scan.c
	$(CC) $(CFLAGS) -c $<

scan.c: scan.l
	$(LEX) -t $< > scan.c

distclean::
	/bin/rm -f parse.c parse.h scan.c parse.output y.output

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(BINOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.o=.d)
endif
