TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

CXXSRC = render.cc oogl.cc object.cc polygons.cc polysphere.cc \
         transform.cc appearance.cc material.cc sphere.cc algebra3.cc \
         color.cc polylines.cc tempinst.cc

CLASSES = Render OOGLRender \
          Transform Appearance Material \
          RenderedObject \
          RenderedSphere RenderedPolygons RenderedPolylines

LIBOBJ = $(CXXSRC:%.cc=%.o)

INC = render.h oogl.h object.h sphere.h color.h \
      algebra3.h parameter.h find.h stack.h polygons.h polysphere.h \
      polylines.h

DEPENDINCLUDE = $(INC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = librender

TESTSRC = rentest.cc
TESTOBJ = $(TESTSRC:%.cc=%.o)
TESTFILES = 

DISTFILES = $(CXXSRC) $(INC) Makefile $(TESTSRC) $(TESTFILES)

default:: $(DEPENDINCLUDE)

interface:: $(DEPENDINCLUDE)

LD = $(CXX)
rentest: $(TESTOBJ) \
	 $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LD) $(LDFLAGS) -o rentest $^ $(SYSLIBS)

rentest.o: rentest.cc
	$(CXX) -DSRCDIR=\"$(SRCDIR)\" $(CPPFLAGS) $(CXXFLAGS) -c $<

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) $(TESTOBJ:%.o=%.d)
endif
