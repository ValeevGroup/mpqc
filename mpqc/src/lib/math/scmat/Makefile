TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
DEFINES += -DSRCDIR=\"$(SRCDIR)\"

default:: $(DEPENDINCLUDE)
default:: libscmat.a

CXXSRC = abstract.cc   local.cc      matrix.cc \
         blkiter.cc    matrix_i.cc   result.cc \
         vector3.cc    vector3_i.cc  block.cc

CSRC = cmatrix.c

LIBOBJ = $(CXXSRC:%.cc=%.o) $(CSRC:%.c=%.o) 

INC =  cmatrix.h matrix_i.h cmatrix.h matrix.h abstract.h blkiter.h local.h \
       result.h vector3.h vector3_i.h block.h

GENINC =

#WEBSRC = blkiter.w  dim.w      include.w  matrix.w   refelem.w \
#         blocks.w   elemop.w   local.w    ref.w


DEPENDINCLUDE = $(INC) $(GENINC)

#$(GENINC): scmat.w $(WEBSRC) 
#	CWEBINPUTS=$(SRCDIR); export CWEBINPUTS; \
#	ctangle $< < /dev/null
#	/bin/rm scmat.c
#	/bin/mv matrix.h.out matrix.h
#	/bin/mv blkiter.h.out blkiter.h
#	/bin/mv local.h.out local.h
#	/bin/mv abstract.h.out abstract.h
#	if cmp matrix.h matrix.h.out > /dev/null; \
#	  then echo header file matrix.h unchanged; \
#	  else (/bin/mv matrix.h.out matrix.h; echo matrix.h updated;) \
#	  fi
#	if cmp blkiter.h blkiter.h.out > /dev/null; \
#	  then echo header file blkiter.h unchanged; \
#	  else (/bin/mv blkiter.h.out blkiter.h; echo blkiter.h updated;) \
#	  fi
#	if cmp local.h local.h.out > /dev/null; \
#	  then echo header file local.h unchanged; \
#	  else (/bin/mv local.h.out local.h; echo local.h updated;) \
#	  fi
#	if cmp abstract.h abstract.h.out > /dev/null; \
#	  then echo header file abstract.h unchanged; \
#	  else (/bin/mv abstract.h.out abstract.h; echo abstract.h updated;) \
#	  fi

#scmat.ps: scmat.dvi
#	dvips -o $@ $^

#see: scmat.dvi
#	xdvi -s 4 -geometry 600x440 $^

#scmat.dvi: scmat.tex
#	tex $^ < /dev/null

#scmat.tex: scmat.w $(WEBSRC) 
#	CWEBINPUTS=$(SRCDIR); export CWEBINPUTS; \
#	cweave $< < /dev/null

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libscmat

TESTSRC = localtest.cc matrixtest.cc
TESTOBJ = matrixtest.o
TESTFILES = # localtest.in

DISTFILES = $(CXXSRC) $(INC) scmat.h $(WEBSRC) Makefile LIBS.h \
            $(TESTSRC) $(TESTFILES)

matrixtest:
	@echo The target must be a specialization such as localtest.

localtest: localtest.o $(TESTOBJ) libscmat.a \
	 $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(CXX) $(CXXFLAGS) -o localtest $^ $(SYSLIBS) 

TAGS: $(CXXSRC) $(CSRC)
	etags $^

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTSRC:%.cc=%.d)
endif

