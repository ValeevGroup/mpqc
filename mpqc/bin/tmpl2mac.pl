#
# This perl script takes a template class definition and converts it to
# a cpp macro.
#

$type = 'notype';
$class = 'noclass';
$macroize = 0;
while (<>) {
    # skip whole line comment and blank lines
    next if (/^\s*(\/\/.*)?$/);
    # strip comments
    s/\/\/.*//;
    if (/^template *< *class +([_A-Za-z][_A-Za-z0-9]*) *>/) { #the template dec
        $type = $1;
    } elsif (/^class +([_A-Za-z][_A-Za-z0-9]*)/) { # the class dec
        $class = $1;
        print "\#define $1_declare($type) \\\n";
        s/^(class +)$class/$1$class \#\# $type/;
        s/\n//;
        # change the names of inherited template classes to macro version
        s/< *$type *>/ \#\# $type/g;
        print "$_ \\\n";
        $macroize = 1;
    } elsif (/^} *; *$/) {      # the end of the class definition
                                # (don't end members like this)
        print "};\n";
        $macroize = 0;
    } elsif (/^ *$/) {          # in case the above misses the definition end
        print "\n";
        $macroize = 0;
    } elsif ($macroize) {       # a line in the body of the class definition
        s/$class *< *$type *>/$class \#\# $type/g;
        s/^( *~? *)$class( *\()/$1$class \#\# $type$2/; # fixes the CTOR's
                                                        # and DTOR's
        # change the names of inherited template classes initializers
        # that appear in CTOR to the macro version
        s/< *$type *>/ \#\# $type/g;
        s/\n//;
        print "$_ \\\n";
    } else {                    # a line outside the class definition
        print "$_";
    }
    
}
