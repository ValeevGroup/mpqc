# -*- Perl -*-

print "\\begin{alltt}\n";
while (<>) {
    if (/-\*- KeyVal -\*-/) {
        # skip the emacs mode line
    }
    else {
        s/<(.*)>/<\\clsnmref{$1}>/;
        print $_;
    }
}
print "\\end{alltt}\n";
