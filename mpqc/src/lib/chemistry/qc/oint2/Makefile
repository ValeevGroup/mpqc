TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

BIN_OR_LIB = LIB
TARGET_TO_MAKE=liboint2

ifeq ($(ARCH),NCUBE_V2)
  CFLAGS = -O3 $(DEFINES) $(INCLUDE)
endif
ifeq ($(ARCH),I860)
  CFLAGS = -O4 -Knoieee -Mnoframe $(DEFINES) $(INCLUDE)
endif


LIBSRC1 = \
int2v0100.c \
int2v0101.c \
int2v0111.c \
int2v0200.c \
int2v0201.c \
int2v0202.c \
int2v0211.c \
int2v0212.c \
int2v0222.c \
int2v1100.c \
int2v1111.c \
int2v1200.c \
int2v1201.c \
int2v1211.c \
int2v1212.c \
int2v1222.c \
int2v2200.c \
int2v2201.c \
int2v2211.c \
int2v2222.c \
int2v0100AB.c \
int2v0101AB.c \
int2v0111AB.c \
int2v0200AB.c \
int2v0201AB.c \
int2v0202AB.c \
int2v0211AB.c \
int2v0212AB.c \
int2v0222AB.c \
int2v1100AB.c \
int2v1111AB.c \
int2v1200AB.c \
int2v1201AB.c \
int2v1211AB.c \
int2v1212AB.c \
int2v1222AB.c \
int2v2200AB.c \
int2v2201AB.c \
int2v2211AB.c \
int2v2222AB.c

LIBSRC2= \
int2v0300.c \
int2v0301.c \
int2v0302.c \
int2v0303.c \
int2v0311.c \
int2v0312.c \
int2v0313.c \
int2v0322.c \
int2v0323.c \
int2v0333.c \
int2v1300.c \
int2v1301.c \
int2v1302.c \
int2v1311.c \
int2v1312.c \
int2v1313.c \
int2v1322.c \
int2v1323.c \
int2v1333.c \
int2v2300.c \
int2v2301.c \
int2v2302.c \
int2v2311.c \
int2v2312.c \
int2v2322.c \
int2v2323.c \
int2v2333.c \
int2v3300.c \
int2v3301.c \
int2v3302.c \
int2v3311.c \
int2v3312.c \
int2v3322.c \
int2v3333.c \
int2v0300AB.c \
int2v0301AB.c \
int2v0302AB.c \
int2v0303AB.c \
int2v0311AB.c \
int2v0312AB.c \
int2v0313AB.c \
int2v0322AB.c \
int2v0323AB.c \
int2v0333AB.c \
int2v1300AB.c \
int2v1301AB.c \
int2v1302AB.c \
int2v1311AB.c \
int2v1312AB.c \
int2v1313AB.c \
int2v1322AB.c \
int2v1323AB.c \
int2v1333AB.c \
int2v2300AB.c \
int2v2301AB.c \
int2v2302AB.c \
int2v2311AB.c \
int2v2312AB.c \
int2v2322AB.c \
int2v2323AB.c \
int2v2333AB.c \
int2v3300AB.c \
int2v3301AB.c \
int2v3302AB.c \
int2v3311AB.c \
int2v3312AB.c \
int2v3322AB.c \
int2v3333AB.c


LIBSRC3=

LIBSRC := $(LIBSRC1)

ifneq ($(ARCH),NCUBE_V2)
  LIBSRC := $(LIBSRC) $(LIBSRC2)

  LIBSRC := $(LIBSRC) $(LIBSRC3)
endif

LIBOBJ = $(LIBSRC:%.c=%.$(OBJSUF))

default:: $(TOPDIR)/lib/liboint2.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/liboint2.so
endif

DISTFILES = Makefile $(LIBSRC1) $(LIBSRC2) $(LIBSRC3)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

# Skip dependency checking for this beast.
#$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
#ifneq ($(DODEPEND),no)
#include $(LIBOBJ:.o=.d)
#endif

