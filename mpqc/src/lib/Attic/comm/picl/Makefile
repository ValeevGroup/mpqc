TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

GENERIC = ext high proc
IPC = msg shm

ifeq ($(DOALLSUBDIRS),yes)
  SUBDIRS = $(GENERIC) $(IPC) paragon n6400 piclnode
else
  SUBDIRS = $(GENERIC)
  ifeq ($(ARCH),SGI)
    SUBDIRS += $(IPC)
  endif
  ifeq ($(ARCH),PARAGON)
    SUBDIRS += paragon
  endif
  ifeq ($(ARCH),NCUBE_V2)
    SUBDIRS += n6400
  endif
  ifeq ($(ARCH),NCUBE_V3)
    SUBDIRS += n6400
  endif
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalSubDirs

