
@libraries = ();
@includes = ();
%defines = ();
%read_files = ();

$debug = 0;

$includes[++$#includes] = ".";
$filename = "";

foreach $arg (@ARGV) {
    if ($arg =~ /^-d$/) {
        $debug = 1;
    }
    elsif ($arg =~ /^-D(.*)$/) {
        local($def) = $1;
        local($symbol) = $1;
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

&process_file($filename);

%current_includes = ();
@libraries = ();
%known_libs = ();
%known_includes = ();
&find_libraries($filename);
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
    local($filename) = shift;
    if ($debug) {
        printf "process_file: filename: %s\n", $filename;
    }

    # find the file
    local($ifile) = "";
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
        print STDERR "listlibs.pl: couldn't find file $file\n";
        exit 1;
    }

    # read the file
    local($filecontents) = "";
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

    # read in other files referenced by this file
    foreach $line (&get_lines($filecontents)) {
        if ($line =~ /^\#\s*include\s*<(.+)>/) {
            local($newfile) = $1;
            if ($read_files{$newfile} eq "") {
                &process_file($newfile);
            }
        }
    }
}

sub get_lines {
    local($filecontents) = shift;
    local(@lines) = ();
    while ($filecontents ne "") {
        $filecontents =~ s/^(.*)\n//;
        $lines[++$#lines] = $1;
    }
    return @lines;
}

sub find_libraries {
    local($filename) = shift;
    if ($current_includes{$filename} == 1) {
        print STDERR "listlibs.pl: recursive include detected for $filename\n";
        exit 1;
    }
    $current_includes{$filename} = 1;
    foreach $line (reverse(get_lines($read_files{$filename}))) {
        if ($line =~ /^\#\s*include\s*<(.+)>/) {
            local($newfile) = $1;
            if ($known_includes{$newfile} != 1) {
                $known_includes{$newfile} = 1;
                &find_libraries($newfile);
            }
        }
        elsif ($known_libs{$line} != 1) {
            $known_libs{$line} = 1;
            $libraries[++$#libraries] = $line;
        }
    }
    delete $current_includes{$filename};
}

sub substitute_defines {
    local($i);
    local($symbol);
    foreach $i (0..$#libraries) {
        foreach $symbol (keys(%defines)) {
            $libraries[$i] =~ s/$symbol/$defines{$symbol}/g;
        }
    }
}

