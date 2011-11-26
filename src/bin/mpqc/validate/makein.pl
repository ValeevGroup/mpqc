#
eval 'exec perl $0 $*'
    if 0;

require QCParse;

$prefix = "";
@files = ();
%files = ();
$writefiles = 1;
$echonames = 0;
$dir = "";
%basissets = (
  "STO-2G" => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],
  "STO-3G" => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36],
  "STO-3G*" => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],
  "STO-6G" => [1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36],
  "MINI (Huzinaga)" => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],
  "MINI (Scaled)" => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],
  "MIDI (Huzinaga)" => [1,2,3,4,5,6,7,8,9,10,11],
  "DZ (Dunning)" => [1,3,5,6,7,8,9,10,13,14,15,16,17],
  "DZP (Dunning)" => [1,3,5,6,7,8,9,10,13,14,15,16,17],
  "DZP + Diffuse (Dunning)" => [1,5,6,7,8,9,10],
  "3-21G" => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36],
  "3-21G*" => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],
  "3-21++G" => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],
  "3-21++G*" => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],
  "4-31G" => [1,2,3,4,5,6,7,8,9,10,15,16,17],
  "6-31G" => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],
  "6-31G*" => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],
  "6-31G**" => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],
  "6-31+G*" => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],
  "6-31++G" => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],
  "6-31++G*" => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],
  "6-31++G**" => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],
  "6-311G" => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,31,32,33,34,35,36],
  "6-311G*" => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,31,32,33,34,35,36],
  "6-311G**" => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,31,32,33,34,35,36],
  "6-311G(2df,2pd)" => [1,2,3,4,5,6,7,8,9,10],
  "6-311++G**" => [1,2,3,4,5,6,7,8,9,10],
  "6-311++G(2d,2p)" => [1,2,3,4,5,6,7,8,9,10],
  "6-311++G(3df,3pd)" => [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],
  "cc-pVDZ" => [1,2,5,6,7,8,9,10,13,14,15,16,17,18],
  "cc-pVTZ" => [1,2,5,6,7,8,9,10,13,14,15,16,17,18],
  "cc-pVQZ" => [1,2,5,6,7,8,9,10,13,14,15,16,17,18],
  "cc-pV5Z" => [1,2,3,4,5,6,7,8,9,10,13,14,15,16,17,18],
  "aug-cc-pVDZ" => [1,2,5,6,7,8,9,10,13,14,15,16,17,18],
  "aug-cc-pVTZ" => [1,2,5,6,7,8,9,10,13,14,15,16,17,18],
  "aug-cc-pVQZ" => [1,2,5,6,7,8,9,10,13,14,15,16,17,18],
  "aug-cc-pV5Z" => [1,2,5,6,7,8,9,10,13,14,15,16,17,18],
  "cc-pCVDZ" => [5,6,7,8,9,10],
  "cc-pCVTZ" => [5,6,7,8,9,10],
  "cc-pCVQZ" => [5,6,7,8,9,10],
  "cc-pCV5Z" => [5,6,7,8,9,10],
  "aug-cc-pCVDZ" => [5,6,7,8,9],
  "aug-cc-pCVTZ" => [5,6,7,8,9,10],
  "aug-cc-pCVQZ" => [5,6,7,8,9,10],
  "aug-cc-pCV5Z" => [5,6,7,8,9],
  "NASA Ames ANO" => [1,5,6,7,8,9,10,13,15,22,26,28],
  "pc-0" => [1,6,7,8,9,14,15,16,17],
  "pc-1" => [1,6,7,8,9,14,15,16,17],
  "pc-2" => [1,6,7,8,9,14,15,16,17],
  "pc-3" => [1,6,7,8,9,14,15,16,17],
  "pc-4" => [1,6,7,8,9,14,15,16,17],
  "pc-0-aug" => [1,6,7,8,9,14,15,16,17],
  "pc-1-aug" => [1,6,7,8,9,14,15,16,17],
  "pc-2-aug" => [1,6,7,8,9,14,15,16,17],
  "pc-3-aug" => [1,6,7,8,9,14,15,16,17],
  "pc-4-aug" => [1,6,7,8,9,14,15,16,17]
);

while ($_ = shift) {
    if (/^-I(.*)/) {
        $prefix = $1;
    }
    elsif (/^-p(.+)$/) {
        $package = $1;
    }
    elsif (/^-p$/) {
        $package = shift;
    }
    elsif (/^-e$/) {
        $writefiles = 0;
        $echonames = 1;
    }
    elsif (/^-d$/) {
        $dir = shift;
        $dir = "$dir/";
    }
    else {
        unshift @files, $_;
    }
}

foreach $i (@files) {
    process_file($i, $package);
}

if ($echonames) {
    for $i (keys(%files)) {
        print "$i\n";
    }
}

sub process_file {
    my $file = shift;
    my $package = shift;
    my $parse = new QCParse;

    my $fullfile = "$prefix/$file";
    $parse->parse_file($fullfile);
    if (! $echonames) {
        print "File: $fullfile\n";
    }

    $file =~ s/^.*\/([^\/]*)$/\1/;
    $file =~ s/\..*$//;

    my $test_vars = {};
    init_var($test_vars, $parse, "basis", "STO-3G");
    init_var($test_vars, $parse, "auxbasis", "");
    init_var($test_vars, $parse, "dfbasis", "");
    init_var($test_vars, $parse, "grid", "default");
    init_var($test_vars, $parse, "symmetry", "C1");
    init_var($test_vars, $parse, "method", "SCF");
    init_var($test_vars, $parse, "calc", "energy");
    init_var($test_vars, $parse, "fzc", 0);
    init_var($test_vars, $parse, "fzv", 0);
    init_var($test_vars, $parse, "docc", "auto");
    init_var($test_vars, $parse, "socc", "auto");
    init_var($test_vars, $parse, "multiplicity", 1);
    init_var($test_vars, $parse, "gradient", "default");
    init_var($test_vars, $parse, "molecule", "molecule");
    init_var($test_vars, $parse, "orthog_method", "default");
    init_var($test_vars, $parse, "lindep_tol", "default");
    init_var($test_vars, $parse, "default_package",
             "MPQC.IntV3EvaluatorFactory");
    # r12 theory params
    init_var($test_vars, $parse, "r12theory", "default");

    my @molecule_symmetry = $parse->value_as_array("test_molecule_symmetry");
    my @molecule_fzc = $parse->value_as_array("test_molecule_fzc");
    my @molecule_fzv = $parse->value_as_array("test_molecule_fzv");
    my @molecule_docc = $parse->value_as_array("test_molecule_docc");
    my @molecule_socc = $parse->value_as_array("test_molecule_socc");
    my @molecule_mult = $parse->value_as_array("test_molecule_multiplicity");
    my @molecule_gradient = $parse->value_as_array("test_molecule_gradient");
    my @molecule_followed = $parse->value_as_array("test_molecule_followed");
    my @molecule_fixed = $parse->value_as_array("test_molecule_fixed");
    my $do_cca = "";
    my $tmp_do_cca = $parse->value("do_cca");
    if ($tmp_do_cca eq "yes") {
         $do_cca = "yes";
    }

    my @keys = keys(%{$test_vars});
    my $index = {};
    my $size = {};
    my $i;
    foreach $i (@keys) {
        $index->{$i} = 0;
        $size->{$i} = scalar(@{$test_vars->{$i}});
    }
    for ($i=0; $i < $nkeys; $i++) { $vals->{$i} = 0; }
    # this is equivalent to a nested loop for each test var,
    # but is much easier to maintain
    do {
        my $basis = $test_vars->{"basis"}->[$index->{"basis"}];
        my $auxbasis = $test_vars->{"auxbasis"}->[$index->{"auxbasis"}];
        my $dfbasis = $test_vars->{"dfbasis"}->[$index->{"dfbasis"}];
        my $grid = $test_vars->{"grid"}->[$index->{"grid"}];
        my $fzc = $test_vars->{"fzc"}->[$index->{"fzc"}];
        my $fzv = $test_vars->{"fzv"}->[$index->{"fzv"}];
        my $docc = $test_vars->{"docc"}->[$index->{"docc"}];
        my $socc = $test_vars->{"socc"}->[$index->{"socc"}];
        my $mult = $test_vars->{"multiplicity"}->[$index->{"multiplicity"}];
        my $gradient = $test_vars->{"gradient"}->[$index->{"gradient"}];
        my $method = $test_vars->{"method"}->[$index->{"method"}];
        my $calc = $test_vars->{"calc"}->[$index->{"calc"}];
        my $symmetry = $test_vars->{"symmetry"}->[$index->{"symmetry"}];
        my $molecule = $test_vars->{"molecule"}->[$index->{"molecule"}];
        my $fixed = $molecule_fixed[$index->{"fixed"}];
        my $followed = $molecule_fixed[$index->{"followed"}];
        my $orthog_method = $test_vars->{"orthog_method"}->[$index->{"orthog_method"}];
        my $lindep_tol = $test_vars->{"lindep_tol"}->[$index->{"lindep_tol"}];
        my $r12theory = $test_vars->{"r12theory"}->[$index->{"r12theory"}];
        my $default_package = $test_vars->{"default_package"}->[$index->{"default_package"}];
        # if i got an array of molecule names then i expect
        # an array of point groups, one for each molecule
        if ($molecule ne "molecule") {
            my $molindex = $index->{"molecule"};
            $symmetry = $molecule_symmetry[$molindex];
            if ($symmetry eq "") {
                printf "\n";
                printf "index = %d\n", $molindex;
                printf "symmetry not set for a molecule array\n";
                exit 1;
            }
            # check for frozen orbitals
            if ($#molecule_fzc >= $molindex) {
                $fzc = $molecule_fzc[$molindex];
            }
            if ($#molecule_fzv >= $molindex) {
                $fzv = $molecule_fzv[$molindex];
            }
            # check for occupations
            if ($#molecule_docc >= $molindex) {
                $docc = $molecule_docc[$molindex];
            }
            if ($#molecule_socc >= $molindex) {
                $socc = $molecule_socc[$molindex];
            }
            if ($#molecule_mult >= $molindex) {
                $mult = $molecule_mult[$molindex];
            }
            if ($#molecule_gradient >= $molindex) {
                $gradient = $molecule_gradient[$molindex];
            }
            # check for fixed coordinates
            $fixed = $molecule_fixed[$molindex];
            if ($fixed eq "-") {
                $fixed = "";
            }
            # check for followed coordinates
            $followed = $molecule_followed[$molindex];
            if ($followed eq "-") {
                $followed = "";
            }
        }

        # the filename to use for the calc
        my $fcalc;

        # only need fzc and fzv for correlated calculations
        if (!($method =~ /MP2/i)
            && !($method =~ /OPT1/i)
            && !($method =~ /OPT2[2]/i)
            && !($method =~ /ZAPT2/i)) {
            $fzc = "";
            $fzv = "";
        }
        if ($calc eq "energy") {
            $parse->set_value("optimize", "no");
            $parse->set_value("frequencies", "no");
            $fcalc = "";
        }
        elsif ($calc eq "opt") {
            $parse->set_value("optimize", "yes");
            $parse->set_value("frequencies", "no");
            $fcalc = "opt";
        }
        elsif ($calc eq "freq") {
            $parse->set_value("optimize", "no");
            $parse->set_value("frequencies", "yes");
            $fcalc = "frq";
        }
        elsif ($calc eq "optfreq") {
            $parse->set_value("optimize", "yes");
            $parse->set_value("frequencies", "yes");
            $fcalc = "optfrq";
        }
        else {
            print "Bad value for calc: $calc\n";
            exit 1;
        }

        my $fextra = ""; # extra filename modifiers
        $parse->set_value("basis", $basis);
        $parse->set_value("auxbasis", $auxbasis);
        $parse->set_value("grid", $grid);
        $parse->set_value("method", $method);
        $parse->set_value("symmetry", $symmetry);
        $parse->set_value("default_package", $default_package);
        $parse->set_value("fzc", $fzc);
        $parse->set_value("fzv", $fzv);
        $parse->set_value("docc", $docc);
        $parse->set_value("socc", $socc);
        $parse->set_value("state", $mult);

	{
	  $dfbasis = "" if ($dfbasis eq "none");
	  $parse->set_value("dfbasis", $dfbasis);
	  $dfbasis = tofilename($dfbasis);
	}

        if ($gradient ne "default") {
            if ($method =~ /v[12](lb)?$/i) {
                # these methods don't support gradients
                $parse->set_value("gradient", "no");
            }
            else {
                $parse->set_value("gradient", $gradient);
            }
        }
        if ($orthog_method ne "default") {
            $parse->set_value("orthog_method", $orthog_method);
            if ($orthog_method eq "gramschmidt") {
                $fextra = "gs$fextra";
            }
            elsif ($orthog_method eq "canonical") {
                $fextra = "can$fextra";
            }
            elsif ($orthog_method eq "symmetric") {
                $fextra = "sym$fextra";
            }
            else {
                $fextra = "$orthog_method$fextra";
            }
        }
        if ($lindep_tol ne "default") {
            $parse->set_value("lindep_tol", $lindep_tol);
            my $ldtolindex = $index->{"lindep_tol"};
            $fextra = "t$ldtolindex$fextra";
        }
	$parse->set_value("r12theory", $r12theory);
        if ($r12theory ne "default") {
            $fextra =  "$r12theory$fextra";
        }
        $parse->set_value("molecule", $parse->value($molecule));
        $parse->set_value("fixed", $parse->value($fixed));
        $parse->set_value("followed", $parse->value($followed));
        my $qcinput = new QCInput($parse);
        my $fmol = $molecule;
        $fmol = "" if ($molecule eq "molecule");
        $fmol = tofilename($fmol);

        # make sure that the basis set exists for all of the
        # atoms in the molecule
        my $molobject = $qcinput->molecule();
        my $allowedatoms = $basissets{$basis};
        my $ok = 1;
        for $symbol (0..($molobject->n_atom()-1)) {
            my $z = $molobject->z($symbol);
            my $gotit = 0;
            for $ztmp (@{$allowedatoms}) {
                if ($ztmp == $z) {
                    $gotit = 1;
                    last;
                }
            }
            $ok = 0 if (! $gotit);
        }

        my $spinok = 1;
        if (($method eq "MP2" || $method eq "mp2" || $method eq "LMP2" || $method eq "lmp2" || $method =~ /MP2V/i) && $mult > 1) {
            $spinok = 0;
        }

        my $inputfile;
        $method = tofilename($method);
        $basis = tofilename($basis);
        $auxbasis = tofilename($auxbasis);
        $symmetry = tofilename($symmetry);
        $fextra = tofilename($fextra);
        if ($do_cca eq "yes"){
             $intpack = tofilename($default_package);
        }
        else {
             $intbuf = "";
             $intpack = "";
        }
        if ($grid eq "default") {$grid = "";}
        my $basename = "$dir$file\_$fmol$method$grid$fzc$fzv$basis$auxbasis$dfbasis$symmetry$fcalc$fextra$intbuf";
        my $writer;

        if ($package eq "g94") {
            $writer = new G94InputWriter($qcinput);
            $inputfile = "$basename.com";
        }
        elsif ($package eq "mpqc") {
            $writer = new MPQCInputWriter($qcinput);
            $inputfile = "$basename.in";
        }

        if (! $ok) {
            if (! $echonames) {
                printf "skipping $inputfile since basis not available\n";
            }
        }
        elsif (! $spinok) {
            if (! $echonames) {
                printf "skipping $inputfile due to mult/method combo\n";
            }
        }
        else {
            if ($writefiles) {
                $writer->write_input("$inputfile.tmp");
                $writer->write_qcinput("$basename.qci");
                if ($files{"$inputfile"}) {
                    unlink("$inputfile.tmp");
                }
                else {
                    $files{"$inputfile"} = 1;
                    my $ret = 1;
                    $ret = system("cmp '$inputfile' '$inputfile.tmp' > /dev/null 2>&1")/256
                        if (-f "$inputfile");
                    if ($ret != 0) {
                        print "writing $inputfile\n";
                        rename("$inputfile.tmp", "$inputfile");
                    }
                    else {
                        unlink("$inputfile.tmp");
                        print "$inputfile is unchanged\n";
                    }
                }
            }
            else {
                $files{"$inputfile"} = 1;
            }
        }

    } while (incr($index,$size));
}

sub incr {
    my $index = shift;
    my $size = shift;
    my @keys = keys(%{$index});
    my $i;
    my $dozero = 0;
    while ($i = shift(@keys)) {
        if ($index->{$i} < $size->{$i} - 1) {
            $index->{$i}++;
            return 1;
        }
        else {
            $index->{$i} = 0;
        }
    }
    return 0;
}

sub init_var {
    my $vars = shift;
    my $parse = shift;
    my $name = shift;
    my $default = shift;
    my $testname = "test_$name";
    my @ar = $parse->value_as_array($testname);
    if ($#ar < 0) { @ar = ( $default ); }
    $vars->{$name} = \@ar;
}

sub tofilename {
    my $raw = shift;
    $raw =~ tr/A-Z/a-z/;
    $raw =~ s/-//g;
    $raw =~ s/ //g;
    $raw =~ s/\*/s/g;
    $raw =~ s/\+/p/g;
    $raw =~ s/\'/prime/g;
    $raw =~ s./.slash.g;
    $raw =~ s/\(/_/g;
    $raw =~ s/\)/_/g;
    $raw;
}
