TOPDIR = ../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libdmt
BIN_OR_LIB = LIB

DEFINES += -DTIMING

# This is for matrix.web.  Adding throttling (on 4-10-92) caused
# some problem with optimization.
ifeq ($(ARCH),NCUBE_V2)
  CFLAGS := $(CFLAGS_DEBUG)
  FFLAGS := $(FFLAGS_DEBUG)
endif

CSRC = hql.c upper.c lower.c ptred.c
FSRC = svdsln.f

INC = libdmt.h matrix.h

DISTFILES = matrix.web $(CSRC) $(FSRC) $(INC) LIBS.h

GENSRC = matrix.c
GENINC = matrix.gbl

LIBCSRC = $(GENSRC) hql.c ptred.c

ifeq ($(MAP_C_F),upper)
  LIBCSRC += upper.c
endif
ifeq ($(MAP_C_F),lower)
  LIBCSRC += lower.c
endif

DEPENDINCLUDE = $(INC) $(GENINC)

LIBOBJ = $(LIBCSRC:%.c=%.$(OBJSUF)) $(FSRC:%.f=%.$(OBJSUF))

###################################################################

default:: $(DEPENDINCLUDE)
default:: libdmt.a

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBCSRC:.c=.d)
endif

##############################################################

ifneq ($(SRCDIR),.)
%.c: %.web
	cp $< .
	$(CTANGLE) $(notdir $<)
	rm -f $(notdir $<)
else
%.c: %.web
	$(CTANGLE) $<
endif

