
--------------------------------------------------------------------------
special variables recognized by makefiles:

DOALLSUBDIRS=yes
  Makes makefiles that descend directories descend all the
directories that they know about, not just those that are
relevant for the particular architecture.  This is useful
for preparing distributions.

DODEPEND=no
  Tells the makefile to not include the .d dependency files.
This is useful for doing clean.

--------------------------------------------------------------------------
variables with special meanings

DISTFILES is the minimal number files needed to rebuild all other files

