TOPDIR = ../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libgeom
BIN_OR_LIB = LIB

#CFLAGS := $(CFLAGS_DEBUG)
#FFLAGS := $(FFLAGS_DEBUG)

WEBSRC = geom.web
GENSRC = $(WEBSRC:%.web=%.c)

INC = libgeom.h retcodes.h
GENINC = geom.gbl

LIBOBJ= $(GENSRC:%.c=%.$(OBJSUF))

DEPENDINCLUDE = $(INC) $(GENINC)

DISTFILES = $(WEBSRC) $(INC) Makefile LIBS.h

###############################################################3

default:: $(DEPENDINCLUDE)
default:: libgeom.a

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

ifneq ($(SRCDIR),.)
%.c: %.web
	cp $< .
	$(CTANGLE) $(notdir $<)
	rm -f $(notdir $<)

%.tex: %.web
	cp $< .
	$(CWEAVE) $(notdir $<)
	rm -f $(notdir $<)
else
%.c: %.web
	$(CTANGLE) $<

%.tex: %.web
	$(CWEAVE) $<
endif
