TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

BIN_OR_LIB = BIN
TARGET_TO_MAKE = x

LIBS = $(LIBDIR)/libsgen.a $(LIBDIR)/libipv2.a

override COPTIONS = -g

SGFILE = test.sgen

SGSRC = \
testinit.c \
testasgn.c \
testip.c \
testiseq.c \
testallc.c \
testzero.c \
testprnt.c

SGINC = $(SGSRC:%.c=%.h)
	
BINOBJ = x.o $(SGSRC:%.c=%.o)

DISTFILES = Makefile x.c $(SGFILE) in.dat LIBS.h out.dat.right

LIBS := $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

default:: x

out.dat: input.dat x
	x $<

.PHONEY: diff
diff: out.dat.right out.dat
	diff $< out.dat

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(SGSRC) $(SGINC) test.h: $(SGFILE)
	../sgen $<

distclean::
	/bin/rm -f $(SGSRC) $(SGINC) test.h

$(BINOBJ:.o=.d): $(SGINC) test.h
ifneq ($(DODEPEND),no)
include $(BINOBJ:.o=.d)
endif
