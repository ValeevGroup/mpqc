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

    local($havedxx) = 0;
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
        elsif ("$file" =~ /\.dxx$/) {
            $havedxx = 1;
            printf DXX "  //\@Include: $dir/$file\n";
        }
    }
}
