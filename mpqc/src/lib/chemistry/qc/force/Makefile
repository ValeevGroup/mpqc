TOPDIR = ../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libforce
BIN_OR_LIB = LIB

DEFINES += -DTIMING

ifeq ($(ARCH),SUN4)
override COPTIONS = -O2
endif

HISTORICALWEBSRC = spcsscf.web mpcsscf.web mposscf.web
WEBSRC = mpcsscf.web mposscf.web
.PRECIOUS: $(WEBSRC:%.web=%.c)
CXXSRC = read.cc

GENSRC = $(WEBSRC:%.web=%.c)
SRCS = $(GENSRC)

INC = libforce.h read.h
GENINC = $(SRCS:%.c=%.gbl)

DEPENDINCLUDE = $(INC) $(GENINC)

LIBOBJ= $(SRCS:%.c=%.$(OBJSUF)) $(CXXSRC:%.cc=%.$(OBJSUF))

################################################################3

default:: $(DEPENDINCLUDE)
default:: $(TOPDIR)/lib/libforce.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libforce.so
endif

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
