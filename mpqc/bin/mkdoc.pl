#!/usr/bin/perl

$srcdir = shift;
$dxx = shift;

if (substr($srcdir,0,1) eq "/") {
    $topdir = "";
}
else {
    $topdir = ".";
}

open(DXX,">$dxx");

&dodir("$srcdir",".",$topdir,0);

close(DXX);


sub dodir {
    local($dir,$objdir,$topdir,$depth) = @_;
    local($file);
    local(@files);

    #print "In directory $dir\n";

    opendir(DIR, $dir) || (warn "Can't open $dir: $!\n", return);
    @files = readdir(DIR);
    closedir(DIR);

    local($havemake) = 0;
    local($haveinc) = 0;
    @includes = ();
    foreach $file (@files) {
        if ($file eq "." || $file eq ".." || $file eq "CVS") {
            next;
        }

        if (-d "$dir/$file") {
            local($nexttop);
            if ($topdir eq ".") {
                $nexttop = "../";
            }
            elsif ($topdir eq "" ) {
                $nexttop = "";
            }
            else {
                $nexttop = "$topdir../";
            }
            &dodir("$dir/$file", "$objdir/$file", $nexttop, $depth+1);
        }
        elsif ("$file" eq Makefile) {
            $havemake = 1;
        }
        elsif ("$file" =~ /\.h$/ && isdocxx("$dir/$file")) {
            $haveinc = 1;
            @includes = (@includes, $file);
        }
    }
    local($incdir) = $haveinc;
    $truncdir = $dir;
    $truncdir =~ s/^\.\.\///;
    if ($incdir) {
        &indent($depth); printf DXX "/**\@name $truncdir\n";
        &indent($depth); printf DXX " */\n";
        &indent($depth); printf DXX "//\@{\n";
    }
    if ($haveinc == 1) {
        foreach $file (sort(@includes)) {
            &indent($depth+1); printf DXX "/**\@name $truncdir/$file\n";
            &indent($depth+1); printf DXX " */\n";
            &indent($depth+1); printf DXX "//\@{\n";
            &indent($depth+2); printf DXX "//\@Include: $dir/$file\n";
            &indent($depth+1); printf DXX "//\@}\n";
        }
    }
    if ($incdir) {
        &indent($depth); printf DXX "//\@}\n";
    }
}

sub isdocxx {
    local($file) = shift;
    open(IN,"<$file");
    while (<IN>) {
        if (/\/\*\*/ || / *\/\/\//) {
            close(IN);
            return 1;
        }
    }
    close(IN);
    return 0;
}

sub indent {
    local($depth) = shift;
    foreach $i (0..$depth) {
        printf DXX "  ";
    }
}
