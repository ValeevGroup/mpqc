# -*- Perl -*-

print "\\begin{alltt}\n";
while (<>) {
    if (/-\*- KeyVal -\*-/) {
        # skip the emacs mode line
    }
    else {
        s/{/\\{/g;
        s/}/\\}/g;
        s/<(.*)>/<\\clsnmref{$1}>/;
        s/\$/\\\$/g;
        print $_;
    }
}
print "\\end{alltt}\n";
