#
eval 'exec perl $0 $*'
    if 0;

require QCResult;

my $file1 = shift;
my $file2 = shift;

check($file1, $file2);

# Takes the name of the output file as the first argument.  It must end in
# a .out. The QCInput file must be in the same directory and must end in a
# .qci.  The optional second argument is the path to an output file that are
# to be compared to the file given by the first argument.
sub check {
    my $fileout = shift;
    my $comparefileout = shift;
    my $file = $fileout;
    $file =~ s/\.out$//;
    my $filein = "$file.qci";
    my $log10 = log(10.0);

    my $result = new QCResult("$filein","$fileout");
    my $ok = "failed";
    if ($result->ok()) {
        if ($result->inputok()) {
            $ok = "ok";
        }
        else {
            $ok = "(ok)";
        }
    }
    else {
        if (! $result->inputok()) {
            $ok = "(failed)";
        }
    }
    $ok = "missing" if (! $result->exists());
    my $basename = $file;
    $basename =~ s=^.*/([^/]*)$=\1=;

    if ($comparefileout eq "") {
        $basename = "$basename:";
        printf "%-28s %s", $basename, $ok;
        if ($result->ok()) {
            printf " %14.8f", $result->energy();
        }
    }
    else {
        my $comparefile = "$comparefileout";
        $comparefile =~ s/\.out$//;
        my $comparebasename = $comparefile;
        $comparebasename =~ s=^.*/([^/]*)$=\1=;
        if ($basename eq $comparebasename) {
            $basename = "$basename:";
            printf "%-28s %s", $basename, $ok;
        }
        else {
            my $files = "$basename/$comparebasename:";
            printf "%-35s %s", $files, $ok;
        }
        if (-f "$comparefile.out") {
            my $comparefileout = "$comparefile.out";
            my $comparefilein = "$comparefile.qci";
            my $cresult = new QCResult($comparefilein,$comparefileout);
            my $compareok = "failed";
            $compareok = "ok" if ($cresult->ok());
            printf " %s", $compareok;
            if ($cresult->ok() && $result->ok()) {
                #printf " %14.8f %14.8f", $result->energy(),$cresult->energy();
                my $diff = abs($result->energy()-$cresult->energy());
                my $ldiff;
                if ($diff == 0) {
                    $ldiff = 99;
                }
                else {
                    $ldiff = -log($diff)/$log10;
                }
                printf " E:%2d", $ldiff;
                print "*" if ($ldiff <= 6);
                if ($result->input()->gradient()
                    && ! $result->input()->optimize()) {
                    my $maxerror = compare_vecs($result->gradient(),
                                                $cresult->gradient());
                    printf " Grad:%2d", $maxerror;
                    print "*" if ($maxerror <= 6);
                }
                if ($result->input()->optimize()) {
                    my $maxerror = compare_vecs(
                                    $result->optmolecule()->geometry(),
                                    $cresult->optmolecule()->geometry());
                    printf " Geom:%2d", $maxerror;
                    print "*" if ($maxerror <= 4);
                }
                if ($result->input()->frequencies()) {
                    my $maxerror = compare_vecs($result->frequencies(),
                                                $cresult->frequencies());
                    printf " Freq:% 2d", $maxerror;
                    print "*" if ($maxerror <= -2);
                }
            }
            else {
                printf " cannot compare since one calc failed";
            }
        }
        else {
            printf " missing";
        }
    }
    printf "\n";
}

sub tofilename {
    my $raw = shift;
    $raw =~ tr/A-Z/a-z/;
    $raw =~ s/-//g;
    $raw =~ s/\*/s/g;
    $raw;
}

sub compare_vecs {
    my $v1ref = shift;
    my $v2ref = shift;
    my @v1 = @{$v1ref};
    my @v2 = @{$v2ref};
    my $e1, $e2;
    my $maxerror = 99;
    my $log10 = log(10.0);
    my $nv1 = @v1;
    my $nv2 = @v2;
    if ($nv1 != $nv2) {
        printf "compare_vecs: vecs not of equal length\n";
        exit 1;
    }
    while (($e1 = shift @v1)
           &&($e2 = shift @v2)) {
        my $diff = abs($e2-$e1);
        my $ldiff;
        if ($diff == 0) {
            $ldiff = 99;
        }
        else {
            $ldiff = -log($diff)/$log10;
        }
        if ($ldiff < $maxerror) { $maxerror = $ldiff; }
    }
    $maxerror;
}

sub print_vec {
    my $v1ref = shift;
    my @v1 = @{$v1ref};
    my $e1;
    while ($e1 = shift @v1) {
        printf " %12.8f\n", $e1;
    }
}
