# Change this to the path to your installed sc-config script.
SCCONFIG = /usr/local/mpqc/current/bin/sc-config
CXX := $(shell $(SCCONFIG) --cxx)
CXXFLAGS := $(shell $(SCCONFIG) --cxxflags)
CPPFLAGS := $(shell $(SCCONFIG) --cppflags)
LIBS := $(shell $(SCCONFIG) --libs)
LIBDIR  := $(shell $(SCCONFIG) --libdir)

mp2: mp2.o
	$(CXX) $(CXXFLAGS) -o $@ $^ -L$(LIBDIR) -lmpqc $(LIBS)
