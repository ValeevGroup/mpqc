#
eval 'exec perl $0 $*'
    if 0;

require QCResult;

my $log10 = log(10.0);

$error = 0;
$refmissing = 0;
$testmissing = 0;
$reffailed = 0;
$testfailed = 0;
$ntest = 0;

if ($ARGV[0] eq "-r") {
    shift;
    $refdir = shift;
    foreach $file1 (@ARGV) {
        $file2 = $file1;
        $file1 =~ s+run+$refdir+;
        check($file1, $file2);
    }
}
elsif ($ARGV[0] eq "-d") {
    shift;
    my $dir = $ARGV[0];
    shift;
    my $rundir = $ARGV[0];
    shift;
    opendir(DIR,"$dir");
    my @files = sort(readdir(DIR));
    closedir(DIR);
    foreach $file (@files) {
        if ($file =~ /.out$/) {
            check("$dir/$file", "$rundir/$file");
        }
    }
}
else {
    my $file1 = shift;
    my $file2 = shift;


# for AIX, which isn't processing the {,} in the argument
    if ($file1 =~ /(.*){(.*),(.*)}(.*)/) {
        $file1 = "$1$2$4";
        $file2 = "$1$3$4";
    }

    check($file1, $file2);
}

print  "*************************************************\n";
printf "* %6d test cases total\n", $ntest;
printf "* %6d numerical discrepancies\n", $error;
printf "* %6d failed reference cases\n", $reffailed;
printf "* %6d missing reference cases\n", $refmissing;
printf "* %6d failed test cases\n", $testfailed;
printf "* %6d missing test cases\n", $testmissing;
print  "*************************************************\n";

if ($error + $testfailed + $testmissing + $reffailed + $refmissing > 0) {
    print "CHECK FAILED\n";
    exit 1;
}
else {
    print "CHECK OK\n";
    exit 0;
}


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
            # use the input file for the reference calculation
            # so it doesn't need to exist in both directories
            my $cresult = new QCResult($filein,$comparefileout);
            my $compareok = "failed";
            $compareok = "ok" if ($cresult->ok());
            printf " %s", $compareok;
            if ($cresult->ok() && $result->ok()) {
                #printf " %14.8f %14.8f", $result->energy(),$cresult->energy();
                my $ldiff = compare_numbers($result->energy(),$cresult->energy());
                printf " E:%2d", $ldiff;
                flagerror() if ($ldiff <= 6);
                if ($result->input()->gradient()
                    && ! $result->input()->optimize()) {
                    my $maxerror = compare_vecs($result->gradient(),
                                                $cresult->gradient());
                    printf " Grad:%2d", $maxerror;
                    flagerror() if ($maxerror <= 6);
                }
                if ($result->input()->optimize()) {
                    my $maxerror = compare_vecs(
                                    $result->optmolecule()->geometry(),
                                    $cresult->optmolecule()->geometry());
                    printf " Geom:%2d", $maxerror;
                    flagerror() if ($maxerror <= 4);
                }
                if ($result->input()->frequencies()) {
                    my $maxerror = compare_vecs($result->frequencies(),
                                                $cresult->frequencies());
                    printf " Freq:% 2d", $maxerror;
                    flagerror() if ($maxerror <= -2);
                }
                if ($result->s2norm() && $cresult->s2norm()) {
                    my $maxerror = compare_numbers($result->s2norm(),
                                                   $cresult->s2norm());
                    printf " S2N:%d", $maxerror;
                    flagerror() if ($maxerror <= 7);
                }
                if (!$cresult->degenerate() &&
                    $result->s2matrix1norm() && $cresult->s2matrix1norm()) {
                    my $maxerror = compare_numbers($result->s2matrix1norm(),
                                                   $cresult->s2matrix1norm());
                    printf " |S2|1:%d", $maxerror;
                    flagerror() if ($maxerror <= 7);
                }
                if ($result->d1mp2() && $cresult->d1mp2()) {
                    my $maxerror = compare_numbers($result->d1mp2(),
                                                   $cresult->d1mp2());
                    printf " D1:%d", $maxerror;
                    flagerror() if ($maxerror <= 7);
                }
                if ($result->d2mp1() && $cresult->d2mp1()) {
                    my $maxerror = compare_numbers($result->d2mp1(),
                                                   $cresult->d2mp1());
                    printf " D2:%d", $maxerror;
                    flagerror() if ($maxerror <= 7);
                }
                if (!$cresult->degenerate() &&
                    $result->s2matrixinfnorm() && $cresult->s2matrixinfnorm()){
                    my $maxerror = compare_numbers($result->s2matrixinfnorm(),
                                                 $cresult->s2matrixinfnorm());
                    printf " |S2|i:%d", $maxerror;
                    flagerror() if ($maxerror <= 7);
                }
                if ($result->npacharge() && $cresult->npacharge()) {
                    my $maxerror = compare_vecs($result->npacharge(),
                                                $cresult->npacharge());
                    printf " NPAq:%d", $maxerror;
                    flagerror() if ($maxerror <= 5);
                    #printf "npacharge\n";
                    #print_vec($result->npacharge());
                }
                if ($result->npashellpop() && $cresult->npashellpop()) {
                    my $maxerror = compare_vecvecs($result->npashellpop(),
                                                   $cresult->npashellpop());
                    printf " NPAp:%d", $maxerror;
                    flagerror() if ($maxerror <= 5);
                    #printf "npashellpop\n";
                    #print_vecvec($result->npashellpop());
                }
                if (!$cresult->degenerate() &&
                    $result->s2large_coef() && $cresult->s2large_coef()) {
                    my $maxerror
                        = compare_vecs_magnitude($result->s2large_coef(),
                                                 $cresult->s2large_coef());
                    printf " S2L:%d", $maxerror;
                    flagerror() if ($maxerror <= 8);
                    my $n = n_nonzero_in_vec($result->s2large_coef());
                    my $xok = compare_string_vecs($result->s2large_i(),
                                                  $cresult->s2large_i(),$n)
                        && compare_string_vecs($result->s2large_a(),
                                               $cresult->s2large_a(),$n);
                    #printf "coef\n";
                    #print_vec($result->s2large_coef());
                    #printf "i\n";
                    #print_string_vec($result->s2large_i());
                    #printf "a\n";
                    #print_string_vec($result->s2large_a());
                    if ($xok) { print " X:OK" }
                    else { print " X:*"; $error++; }
                }
                if (!$cresult->degenerate() &&
                    $result->d1large_coef() && $cresult->d1large_coef()) {
                    my $maxerror
                        = compare_vecs_magnitude($result->d1large_coef(),
                                                 $cresult->d1large_coef());
                    printf " D1L:%d", $maxerror;
                    flagerror() if ($maxerror <= 7);
                    my $n = n_nonzero_in_vec($result->d1large_coef());
                    my $xok = compare_string_vecs($result->d1large_i(),
                                                  $cresult->d1large_i(),$n)
                        && compare_string_vecs($result->d1large_j(),
                                               $cresult->d1large_j(),$n)
                        && compare_string_vecs($result->d1large_a(),
                                               $cresult->d1large_a(),$n)
                        && compare_string_vecs($result->d1large_b(),
                                               $cresult->d1large_b(),$n)
                        && compare_string_vecs($result->d1large_spin(),
                                               $cresult->d1large_spin(),$n);
                    if ($xok) { print " X:OK" }
                    else { print " X:*"; $error++; }
                    #printf "coef\n";
                    #print_vec($result->d1large_coef());
                    #printf "i\n";
                    #print_string_vec($result->d1large_i());
                    #printf "j\n";
                    #print_string_vec($result->d1large_j());
                    #printf "a\n";
                    #print_string_vec($result->d1large_a());
                    #printf "b\n";
                    #print_string_vec($result->d1large_b());
                    #printf "spin\n";
                    #print_string_vec($result->d1large_spin());
                }
            }
            else {
                if (($result->exists() && $cresult->exists())
                    ||($result->exists() && !$result->ok())
                    ||($cresult->exists() && !$cresult->ok())) {
                    printf " cannot compare since one calc failed";
                }
                if (!$result->exists()) {
                    $refmissing++;
                }
                elsif (!$result->ok()) {
                    $reffailed++;
                }
                if (!$cresult->exists()) {
                    $testmissing++;
                }
                elsif (!$cresult->ok()) {
                    $testfailed++;
                }
            }
        }
        else {
            printf " missing";
            $testmissing++;
        }
    }
    $ntest++;
    printf "\n";
}

sub flagerror {
    print "*";
    $error++;
}

sub tofilename {
    my $raw = shift;
    $raw =~ tr/A-Z/a-z/;
    $raw =~ s/-//g;
    $raw =~ s/\*/s/g;
    $raw;
}

sub compare_numbers {
    my $num1 = shift;
    my $num2 = shift;
    my $diff = abs($num1-$num2);
    my $ldiff;
    if ($diff == 0) {
        $ldiff = 99;
    }
    else {
        $ldiff = -log($diff)/$log10;
    }
    $ldiff;
}

# counts how many elements until we get to the first
# element equal to zero
sub n_nonzero_in_vec {
    my $vref = shift;
    my @v = @{$vref};
    my $n = 0;
    my $e1;
    while (($e1 = shift @v1)) {
        last if (abs($e1) < 1.0e-6);
        $n = $n + 1;
    }
    $n;
}

sub compare_vecs {
    my $v1ref = shift;
    my $v2ref = shift;
    my @v1 = @{$v1ref};
    my @v2 = @{$v2ref};
    my $e1, $e2;
    my $maxerror = 99;
    my $nv1 = @v1;
    my $nv2 = @v2;
    if ($nv1 != $nv2) {
        printf "<compare_vecs: vecs not of equal length>";
        return -$maxerror;
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

sub compare_vecs_magnitude {
    my $v1ref = shift;
    my $v2ref = shift;
    my @v1 = @{$v1ref};
    my @v2 = @{$v2ref};
    my $e1, $e2;
    my $maxerror = 99;
    my $nv1 = @v1;
    my $nv2 = @v2;
    if ($nv1 != $nv2) {
        printf "<compare_vecs_magnitude: vecs not of equal length>";
        return -$maxerror;
    }
    while (($e1 = shift @v1)
           &&($e2 = shift @v2)) {
        my $diff = abs(abs($e2)-abs($e1));
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

sub compare_vecvecs {
    my $v1ref = shift;
    my $v2ref = shift;
    my @v1 = @{$v1ref};
    my @v2 = @{$v2ref};
    my $e1, $e2;
    my $maxerror = 99;
    my $nv1 = @v1;
    my $nv2 = @v2;
    if ($nv1 != $nv2) {
        printf "<compare_vecvecs: vecs not of equal length>";
        return -$maxerror;
    }
    while (($e1 = shift @v1)
           &&($e2 = shift @v2)) {
        my $diff = abs($e2-$e1);
        my $ldiff = compare_vecs($e1,$e2);
        if ($ldiff < $maxerror) { $maxerror = $ldiff; }
    }
    $maxerror;
}

# returns 1 if the vecs are identical for as many elements
# are given in the third argument
sub compare_string_vecs {
    my $v1ref = shift;
    my $v2ref = shift;
    my $n = shift;
    my @v1 = @{$v1ref};
    my @v2 = @{$v2ref};
    my $nv1 = @v1;
    my $nv2 = @v2;
    if ($nv1 != $nv2) {
        printf "<compare_vecs: vecs not of equal length>";
        return 0;
    }
    my $e1, $e2;
    my $i = 0;
    while (($e1 = shift @v1)
           &&($e2 = shift @v2) && $i < $n) {
        if ($e1 ne $e2) { return 0; }
        $i = $i + 1;
    }
    1;
}

sub print_vec {
    my $v1ref = shift;
    my @v1 = @{$v1ref};
    my $e1;
    while ($e1 = shift @v1) {
        printf " %12.8f\n", $e1;
    }
}

sub print_vecvec {
    my $v1ref = shift;
    my @v1 = @{$v1ref};
    my $e1;
    while ($e1 = shift @v1) {
        my @v2 = @{$e1};
        my $e2;
        while ($e2 = shift @v2) {
            printf " %12.8f", $e2;
        }
        printf "\n";
    }
}

sub print_string_vec {
    my $v1ref = shift;
    my @v1 = @{$v1ref};
    my $e1;
    while ($e1 = shift @v1) {
        printf " %s\n", $e1;
    }
}
