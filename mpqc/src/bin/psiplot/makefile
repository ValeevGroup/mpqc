TOPDIR = ../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = psiplot
#BIN_OR_LIB = BIN

#CC = gcc
#CINCLUDE = -I/usr/explorer/include -I/net/munin/usr/people/cljanss/include
#C++INCLUDE = $(CINCLUDE)
#CFLAGS = -g $(CINCLUDE)
#CXX = g++
#CXXFLAGS = -g $(C++INCLUDE)

INCLUDE := $(INCLUDE) -I/usr/explorer/include
CXXINCLUDE := $(CXXINCLUDE) -I/usr/explorer/include

CXGENMDW = /usr/explorer/lib/cxMdw

MOD_MAIN_OBJ = /usr/explorer/lib/cxMcw.o

default:: psiplot

BINOBJ = psiplot.cx.o psiplot.mdw.o psiplot.o

EXLIBS = -L/lib -L/usr/lib  \
         -L/usr/explorer/lib -L.  -L/usr/explorer/lib -lfortranApi -lapi \
         -lcxna -L/usr/lib -lm   -lsun -lmalloc -lmpc -lc_s

psiplot:  $(BINOBJ) \
	 $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(CXX) $(CXXFLAGS) -o $@  ${MOD_MAIN_OBJ} $^ \
	$(EXLIBS)
#	 -L/usr/explorer/lib -lapi -lcx_s \
#	 -lm -lsun -lc_s -lmpc

psiplot.mdw.c:	 psiplot.mres ${CXGENMDW}
	/bin/cp $< /tmp/cxmdw_tmp.mres
	(cd /tmp; ${CXGENMDW} -mdw cxmdw_tmp)
	/bin/cp /tmp/cxmdw_tmp.mdw.c $@
	/bin/rm /tmp/cxmdw_tmp.*
#	${CXGENMDW} -d $(<D) $(basename $(<F))

psiplot.cx.c:	 psiplot.mres ${CXGENMDW}
	/bin/cp $< /tmp/cxcx_tmp.mres
	(cd /tmp; ${CXGENMDW} -cx cxcx_tmp)
	/bin/cp /tmp/cxcx_tmp.cx.c $@
	/bin/rm /tmp/cxcx_tmp.*

EXPLORERDIR = /net/munin/usr/people/cljanss/explorer/modules
install:: install_explorer install_mres
install_explorer: psiplot
	cp $^ $(EXPLORERDIR)
install_mres: psiplot.mres
	cp $^ $(EXPLORERDIR)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(BINOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.o=.d)
endif
