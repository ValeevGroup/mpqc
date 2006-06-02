
@libraries = ();
@includes = ();
%defines = ();
%read_files = ();

$debug = 0;

$includes[++$#includes] = ".";
$filename = "";

@processing = ("Processing files:");

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
        $includes[++$#includes] = $1;
    }
    else {
        $filename = $arg;
    }
}

if ($filename eq "") {
    print STDERR "listlibs.pl: require a filename\n";
    exit 1;
}

if (-d "$filename") {
    %current_includes = ();
    @libraries = ();
    %known_libs = ();
    %known_includes = ();
    &process_directory($filename);
}
else {
    &process_file($filename);
    %current_includes = ();
    @libraries = ();
    %known_libs = ();
    %known_includes = ();
    &find_libraries($filename);
}

@libraries = reverse(@libraries);

print "got $#libraries of them\n" if ($debug);

&substitute_defines();

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
        foreach my $include (@includes) {
            $ifile = "$include/$filename";
            if ($debug) {
                #printf STDERR "process_file: looking for: %s\n", $ifile;
            }
            if (-f $ifile) { last; }
        }
    }
    if ($ifile eq "" || ! -f $ifile) {
        print STDERR "listlibs.pl: couldn't find file \"$ifile\"\n";
        for (my $i=0; $i <= $#processing; $i++) {
          printf STDERR "  %s\n", "$processing[$i]";
        }
        exit 1;
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
    # an empty file will look like a new file below so put in a newline
    if ($read_files{$filename} eq "") {
        $read_files{$filename} = "\n"
    }

    # process the lines of the file
    $processing[++$#processing] = "$ifile";
    get_lines($filecontents);
    --$#processing;
}

sub get_lines {
    my $filecontents = shift;
    my @lines = ();
    my $ifdepth = 0;
    while ($filecontents ne "") {
        # get next line
        $filecontents =~ s/^(.*)\n//;
        my $line = $1;
        # remove comments
        $line =~ s/\/\/.*$//;
        # remove leading trailing whitespace
        $line =~ s/^\s*//;
        $line =~ s/\s*$//;
        # this only handles ifdef's that are one level deep
        if ($line =~ /\#\s*(ifn*def)\s+([a-zA-Z_]\w*)/) {
            my $command = $1;
            my $symbol = $2;
            if (defined $defines{$symbol}) {
                print STDERR "\$defines{$symbol} = $defines{$symbol}\n"
                    if ($debug);
            }
            else {
                print STDERR "$symbol is not defined\n"
                    if ($debug);
            }
            if (($command eq "ifdef" && ! defined $defines{$symbol})
                ||($command eq "ifndef" && defined $defines{$symbol})) {
                while ($filecontents ne "") {
                    $filecontents =~ s/^(.*)\n//;
                    my $tline = $1;
                    last if ($tline =~ /\#\s*endif/);
                }
            }
        }
        elsif ($line =~ /\#\s*endif/) {
            # eat this endif
        }
        elsif ($line =~ /^\#\s*include\s*<(.+)>/) {
            my $newfile = $1;
            # keep this line for find_libraries
            $lines[++$#lines] = $line;
            if ($read_files{$newfile} eq "") {
                &process_file($newfile);
            }
        }
        elsif ($line =~ /^\#\s*define\s+(\w+)\s+(\w+)/) {
            my $sym = $1;
            my $value = $2;
            $defines{$sym} = $value;
            print STDERR "\$defines{$sym} = $defines{$sym}\n" if ($debug);
        }
        elsif ($line =~ /^\#\s*define\s+(\w+)/) {
            my $sym = $1;
            my $value = 1;
            $defines{$sym} = $value;
            print STDERR "\$defines{$sym} = $defines{$sym}\n" if ($debug);
        }
        else {
            if ($line ne "") {
                $lines[++$#lines] = $line;
                print STDERR "got line: $line\n" if ($debug);
            }
        }
    }
    return @lines;
}

sub find_libraries {
    my $filename = shift;
    if ($debug) {
        printf STDERR "find_libraries: filename: $filename\n";
    }
    if ($current_includes{$filename} == 1) {
        print STDERR "listlibs.pl: recursive include detected for $filename\n";
        exit 1;
    }
    $current_includes{$filename} = 1;
    foreach my $line (reverse(&get_lines($read_files{$filename}))) {
        if ($debug) {
            printf STDERR "find_libraries: line: $line\n";
        }
        if ($line =~ /^\#\s*include\s*<(.+)>/) {
            my $newfile = $1;
            if ($known_includes{$newfile} != 1) {
                $known_includes{$newfile} = 1;
                &find_libraries($newfile);
            }
        }
        elsif ($line =~ /^\#\s*define\s+/) {
            # skip this line
        }
        elsif ($known_libs{$line} != 1) {
            $known_libs{$line} = 1;
            $libraries[++$#libraries] = $line;
        }
    }
    delete $current_includes{$filename};
}

sub substitute_defines {
    foreach my $i (0..$#libraries) {
        foreach my $symbol (keys(%defines)) {
            $libraries[$i] =~ s/$symbol/$defines{$symbol}/g;
        }
    }
}

sub process_directory {
    my $dir = shift;
    opendir(DIR,"$dir");
    my @files = readdir(DIR);
    closedir(DIR);
    foreach my $i (@files) {
        if ("$i" eq "." || "$i" eq "..") {
            # skip
        }
        elsif (-d "$dir/$i") {
            process_directory("$dir/$i");
        }
        elsif ("$i" eq "LIBS.h") {
            process_file("$dir/$i");
            &find_libraries("$dir/$i");
        }
    }
}
