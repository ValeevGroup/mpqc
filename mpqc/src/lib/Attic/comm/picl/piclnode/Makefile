TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libpiclnode
BIN_OR_LIB = LIB

INCLUDE += -I$(SRCDIR)

MACHD =

ifeq ($(ARCH),I860)
  MACHD = i860
  INCLUDE += -I$(SRCDIR) -I$(SRCDIR)/i860
endif


PICLSRC = \
  antipode0.c barrier0.c bcast0.c bcast1.c check0.c \
  clocksync0.c close0.c datasize0.c envclocksync.c enverrmsg.c \
  envflush.c envnflush.c gand0.c gather0.c gcomb0.c getarc0.c ginv0.c \
  gmax0.c gmin0.c gor0.c gprod0.c gray0.c gsum0.c gxor0.c \
  host0.c malloc0.c message0.c neighbor0.c open0.c probe0.c \
  recv0.c recvbegin0.c recvend0.c recvinfo0.c send0.c sendbegin0.c \
  sendend0.c setarc0.c sync0.c traceblkbeg.c traceblkend.c \
  traceblock.c traceexit.c tracefiles.c traceflush.c traceinfo.c \
  tracelevel.c tracemark.c tracemarks.c tracemsg.c tracenode.c who0.c

ifdef MACHD
MACHSRC = $(MACHD)/clock0.c $(MACHD)/envhost.c $(MACHD)/envnodeidle.c \
          $(MACHD)/envnodetime.c $(MACHD)/envopen.c $(MACHD)/envprobe.c \
          $(MACHD)/envrecv.c $(MACHD)/envrecvbeg.c $(MACHD)/envrecvend.c \
          $(MACHD)/envsend.c $(MACHD)/envsendbeg.c $(MACHD)/envsendend.c \
          $(MACHD)/envsync.c $(MACHD)/f2c.c
endif

PICLOBJ = $(PICLSRC:%.c=%.$(OBJSUF))
MACHOBJ = $(MACHSRC:$(MACHD)/%.c=%.$(OBJSUF))

LIBSRC = $(PICLSRC) $(MACHSRC)
LIBOBJ = $(PICLOBJ) $(MACHOBJ)

###############################################################

default:: libpiclnode.a

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif


%.o: $(SRCDIR)/$(MACHD)/%.c
	$(CC) $(CFLAGS) $< -c $(OUTPUT_OPTION)

%.d: $(SRCDIR)/$(MACHD)/%.c
	$(CC) -M $(CFLAGS) $< | sed 's/$*.o/$*.o $*.d/g' > $(@F)

