
@libraries = ();
@includes = ();
%defines = ();
%read_files = ();

$debug = 0;

push @includes, ".";
$filename = "";

foreach $arg (@ARGV) {
    if ($arg =~ /^-d$/) {
        $debug = 1;
    }
    elsif ($arg =~ /^-D(.*)$/) {
        my $def = $1;
        my $symbol = $1;
        $def =~ s/^.*=//;
        $symbol =~ s/=.*$//;
        $defines{$symbol} = $def;
    }
    elsif ($arg =~ /^-I(.*)$/) {
        push @includes, $1;
    }
    else {
        $filename = $arg;
    }
}

if ($filename eq "") {
    die "listlibs requires a filename";
}

process_file($filename);

%current_includes = ();
@libraries = ();
%known_libs = ();
%known_includes = ();
find_libraries($filename);

print "got $#libraries of them\n" if ($debug);

substitute_defines();

foreach $i (0..$#libraries) {
    printf "%s", $libraries[$i];
    if ($i < $#libraries) { printf " "; }
}
printf "\n";

###########################################################################

sub process_file {
    my $filename = shift;
    if ($debug) {
        printf "process_file: filename: %s\n", $filename;
    }

    # find the file
    my $ifile = "";
    if ($filename =~ /^\//) {
        $ifile = $filename;
    }
    else {
        foreach $include (@includes) {
            $ifile = "$include/$filename";
            if ($debug) {
                #printf "process_file: looking for: %s\n", $ifile;
            }
            if (-f $ifile) { last; }
        }
    }
    if ($ifile eq "" || ! -f $ifile) {
        die "couldn't find file $file";
    }

    # read the file
    my $filecontents = "";
    open(IFILE,"<$ifile");
    while (<IFILE>) {
        if (/^\s*$/) { next; }
        $filecontents = "$filecontents$_";
    }
    close(IFILE);
    $read_files{$filename} = $filecontents;

    # read in other files referenced by this file
    foreach $line (get_lines($filecontents)) {
        if ($line =~ /^\#\s*include\s*<(.+)>/) {
            my $newfile = $1;
            if (!exists($read_files{$newfile})) {
                process_file($newfile);
            }
        }
    }
}

sub get_lines {
    my $filecontents = shift;
    my @lines = ();
    while ($filecontents ne "") {
        $filecontents =~ s/^(.*)\n//;
        push @lines, $1;
    }
    return @lines;
}

sub find_libraries {
    my $filename = shift;
    if (exists($current_includes{$filename})) {
        die "recursive included detected for $filename";
    }
    $current_includes{$filename} = 1;
    foreach $line (reverse(get_lines($read_files{$filename}))) {
        if ($line =~ /^\#\s*include\s*<(.+)>/) {
            my $newfile = $1;
            if (!exists($known_includes{$newfile})) {
                $known_includes{$newfile} = 1;
                find_libraries($newfile);
            }
        }
        elsif (!exists($known_libs{$line})) {
            $known_libs{$line} = 1;
            unshift @libraries, $line;
        }
    }
    delete $current_includes{$filename};
}

sub substitute_defines {
    my $i;
    my $symbol;
    foreach $i (0..$#libraries) {
        foreach $symbol (keys(%defines)) {
            $libraries[$i] =~ s/$symbol/$defines{$symbol}/g;
        }
    }
}

