#
eval 'exec perl $0 $*'
    if 0;

require QCParse;
require Molecule;

##########################################################################

package QCResult;

# this seems to not work as expected
$fltrx = "([-+]?(?:\\d+\\.\\d*|\\d+|\\.\\d+)(?:[eEdD][+-]?\\d+)?)";

$have_nodenum = 0;

sub test {
    my $i;
    @nums = ("-1.0", "-1", "-.1", "-1.",
             "-1.0E-04", "-1E-04", "-.1E-04", "-1.E-04");
    foreach $i ( @nums ) {
        $i =~ /$fltrx/;
        my $flt = $1;
        printf "%10s %10s %12.8f\n", $i, $flt, $flt;
    }
}

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {};
    bless $self, $class;
    $self->initialize(@_);
    return $self;
}

sub initialize {
    my $self = shift;
    my $infile = shift;
    my $outfile = shift;
    my $parse = new QCParse();
    $parse->parse_file($infile);
    $self->{"qcinput"} = QCInput::new QCInput($parse);

    $self->{"exists"} = 0;
    $self->{"ok"} = 0;
    if (-e $outfile) {
        $self->{"exists"} = 1;
        open(OUTFILE,"<$outfile");
        my $first = 1;
        while (<OUTFILE>) {
            if ($first && /^ *[0-9]+:/) {
                $have_nodenum = 1;
                $first = 0;
            }
            s/^ *[0-9]+:// if ($have_nodenum);
            if (/^\s*MPQC:/) {
                $self->parse_mpqc(\*OUTFILE);
                last;
            }
            elsif (/^\s*Entering Gaussian System/) {
                $self->parse_g94(\*OUTFILE);
                last;
            }
        }
        close(OUTFILE);
    }
}

sub parse_g94 {
    my $self = shift;
    my $out = shift;
    my $scfenergy = "";
    my $mp2energy = "";
    my $ccsd_tenergy = "";
    my $t1norm = "";
    my $optconverged = 0;
    my $molecule;
    my $havefreq = 0;
    my $freq = [];
    my $ifreq = 0;
    my $b3pw91energy = "";
    while (<$out>) {
        s/^ *[0-9]+:// if ($have_nodenum);
        if (/^\s*SCF Done:  E\(RHF\) =\s*$fltrx\s/) {
            $scfenergy = $1;
        }
        elsif (/^\s*E2\s*=\s*$fltrx\s*EUMP2\s*=\s*$fltrx/) {
            $mp2energy = $2;
        }
        elsif (/^\s*CCSD\(T\)\s*=\s*$fltrx/) {
            $ccsd_tenergy = $1;
            $ccsd_tenergy =~ s/[DdE]/e/;
        }
        elsif (/E\(RB\+HF-PW91\)\s*=\s*$fltrx/) {
            $b3pw91energy = $1;
            $b3pw91energy =~ s/[DdE]/e/;
        }
        elsif (/^\s*T1 Diagnostic\s*=\s*$fltrx/) {
            $t1norm = $1;
        }
        elsif (/^\s*-- Stationary point found./
               || /CONVERGENCE CRITERIA APPARENTLY SATISFIED/) {
            $optconverged = 1;
        }
        elsif ($optconverged
               && (/^\s*Input orientation:/ || /^\s*Standard orientation:/)) {
            <$out>; <$out>; <$out>; <$out>;
            my $molstr = "";
            while (<$out>) {
                if (! /^\s+\d+\s+(\d+)\s+$fltrx\s+$fltrx\s+$fltrx\s+/) {
                    last;
                }
                $molstr = "${molstr} $1 $2 $3 $4\n";
            }
            $molecule = new Molecule($molstr);
        }
        elsif (/^\s*Frequencies --\s+$fltrx\s+$fltrx\s+$fltrx/) {
            $freq->[$ifreq] = $1; $ifreq++;
            $freq->[$ifreq] = $2; $ifreq++;
            $freq->[$ifreq] = $3; $ifreq++;
            $havefreq = 1;
        }
        elsif (/^\s*Frequencies --\s+$fltrx\s+$fltrx/) {
            $freq->[$ifreq] = $1; $ifreq++;
            $freq->[$ifreq] = $2; $ifreq++;
            $havefreq = 1;
        }
        elsif (/^\s*Frequencies --\s+$fltrx\s/) {
            $freq->[$ifreq] = $1; $ifreq++;
            $havefreq = 1;
        }
    }
    $self->{"scfenergy"} = $scfenergy;
    $self->{"mp2energy"} = $mp2energy;
    $self->{"b3pw91energy"} = $b3pw91energy;
    $self->{"optconverged"} = $optconverged;
    if ($optconverged) {
        $self->{"optmolecule"} = $molecule;
    }
    $self->{"have_frequencies"} = $havefreq;
    if ($havefreq) {
        my @tmp = sort(@$freq);
        $freq = \@tmp;
        $self->{"freq"} = $freq;
    }
    my $qcinput = $self->{"qcinput"};
    my $method = $qcinput->method();
    if ($method eq "MP2") {
        $self->{"energy"} = $mp2energy;
    }
    elsif ($method eq "CCSD(T)") {
        $self->{"energy"} = $ccsd_tenergy;
    }
    elsif ($method eq "B3-PW91") {
        $self->{"energy"} = $b3pw91energy;
    }
    elsif ($method eq "SCF"
           || $method eq "ROSCF") {
        $self->{"energy"} = $scfenergy;
    }

    $self->{"t1norm"} = $t1norm;

    $self->{"ok"} = 0;
    if ($self->{"energy"} ne "") {
        if ($qcinput->optimize()) {
            if ($self->{"optconverged"}) {
                $self->{"ok"} = 1
                }
            else {
                #printf "not ok because not converged\n";
            }
        }
        else {
            $self->{"ok"} = 1;
        }
        if ($qcinput->frequencies() && ! $havefreq) {
            $self->{"ok"} = 0;
            #printf "not ok because no freq\n";
        }
    }
    else {
        #printf "not ok because no energy\n";
    }
}

sub parse_mpqc {
    my $self = shift;
    my $out = shift;
    my $scfenergy = "";
    my $opt1energy = "";
    my $opt2energy = "";
    my $zapt2energy = "";
    my $mp2energy = "";
    my $optconverged = 0;
    my $molecule;
    my $havefreq = 0;
    my $freq;
    my $wante = 1;
    my $ifreq = 0;
    my $grad;
    my $ngrad;
    my $state = "none";
    my $molecularenergy = "";
    my $s2norm = "";
    my $to_angstrom = 0.52917706;
    my $geometry_conversion = $to_angstrom;
    my $error = "";
    my $ccsd_energy = "";
    my $ccsd_t_energy = "";
    my $t1norm = "";
    my $t12norm = "";
    my $psioutput = "";
    while (<$out>) {
        s/^ *[0-9]+:// if ($have_nodenum);
        if ($state eq "read grad" && $wante) {
            if (/^\s*([0-9]+\s+[A-Za-z]+\s+)?([\-\.0-9]+)\s+([\-\.0-9]+)\s+([\-\.0-9]+)/) {
                $grad->[$ngrad + 0] = $2;
                $grad->[$ngrad + 1] = $3;
                $grad->[$ngrad + 2] = $4;
                $ngrad = $ngrad + 3;
            }
            elsif (/^\s*([0-9]+\s+)?([\-\.0-9]+)/) {
                $grad->[$ngrad] = $2;
                $ngrad = $ngrad + 1;
            }
            else {
                $self->{"grad"} = $grad;
                $state = "none";
            }
        }

        if ($wante && /total scf energy =\s+$fltrx/) {
            $scfenergy = $1;
        }
        elsif ($wante && /OPT1 energy .*:\s+$fltrx/) {
            $opt1energy = $1;
        }
        elsif ($wante && /OPT2 energy .*:\s+$fltrx/) {
            $opt2energy = $1;
        }
        elsif ($wante && /ZAPT2 energy .*:\s+$fltrx/) {
            $zapt2energy = $1;
        }
        elsif ($wante && /MP2 energy .*:\s+$fltrx/) {
            $mp2energy = $1;
        }
        elsif (/CSCF: An SCF program written in C/) {
            $psioutput = 1;
        }
        elsif ($psioutput && $wante
               && /total energy       =\s*$fltrx/) {
            $scfenergy = $1;
        }
        elsif ($psioutput && $wante
               && /\s*ITER\s+CORRELATION ENERGY\s+T1 2-NORM\s+T1 DIAG/) {
            # this is a PSI CC output embedded in an MPQC output
            # grab iteration 0
            my $ecorr;
            <$out>;
            while (<$out>) {
                if (/^\s*\d+\s+$fltrx\s+$fltrx\s+$fltrx/) {
                    $ecorr = $1;
                    $t12norm = $2;
                    $t1norm = $3;
                }
                else {
                    last;
                }
            }
            $ccsd_energy = $ecorr + $scfenergy;
        }
# IMBN: Added the following elsif to handle new cc output 
        elsif ($psioutput && $wante
               && /\s*ITER\s+CORRELATION ENERGY\s+T1 DIAG\s+D1\(CCSD\)/) {
            # this is a PSI CC output embedded in an MPQC output
            # grab iteration 0
            my $ecorr;
            <$out>;
            while (<$out>) {
                if (/^\s*\d+\s+$fltrx\s+$fltrx\s+$fltrx/) {
                    $ecorr = $1;
                    $t1norm = $2;
                    $t12norm = $3;
                }
                else {
                    last;
                }
            }
            $ccsd_energy = $ecorr + $scfenergy;
        }
        elsif ($wante && /Value of the MolecularEnergy:\s+$fltrx/) {
            $molecularenergy = $1;
        }
        elsif ($wante && /S2\(ov\) matrix 1-norm\s*=\s*$fltrx/) {
            $self->{"s2matrix1norm"} = $1;
        }
        elsif ($wante && /There are degenerate orbitals within an irrep/) {
            $self->{"degenerate"} = 1;
        }
        elsif ($wante && /D1\(MP2\)\s*=\s*$fltrx/) {
            $self->{"d1mp2"} = $1;
        }
        elsif ($wante && /D2\(MP1\)\s*=\s*$fltrx/) {
            $self->{"d2mp1"} = $1;
        }
        elsif ($wante && /S2\(ov\) matrix inf-norm\s*=\s*$fltrx/) {
            $self->{"s2matrixinfnorm"} = $1;
        }
        elsif ($wante && /S2 norm\s*=\s*$fltrx/) {
            $s2norm = $1;
        }
        elsif ($wante && /Largest S2 values.*:/) {
            my $s2large_coef = [];
            my $s2large_i = [];
            my $s2large_a = [];
            my $is2 = 0;
            while (<$out>) {
                if (/^\s*([0-9]+)\s+$fltrx\s+([0-9]+)\s+([A-Za-z0-9\'\"]+)\s+->\s+([0-9]+)\s+([A-Za-z0-9\'\"]+)\s*$/) {
                    $s2large_coef->[$is2] = $2;
                    $s2large_i->[$is2] = "$3$4";
                    $s2large_a->[$is2] = "$5$6";
                    $is2 = $is2 + 1;
                }
                else {
                    last;
                }
            }
            $self->{"s2large_coef"} = $s2large_coef;
            $self->{"s2large_i"} = $s2large_i;
            $self->{"s2large_a"} = $s2large_a;
        }
        elsif ($wante && /Largest first order coefficients.*:/) {
            my $d1large_coef = [];
            my $d1large_i = [];
            my $d1large_j = [];
            my $d1large_a = [];
            my $d1large_b = [];
            my $d1large_spin = [];
            my $id1 = 0;
            while (<$out>) {
                if (/^\s*([0-9]+)\s+$fltrx\s+([0-9]+)\s+([A-Za-z0-9\'\"]+)\s+([0-9]+)\s+([A-Za-z0-9\'\"]+)\s+->\s+([0-9]+)\s+([A-Za-z0-9\'\"]+)\s+([0-9]+)\s+([A-Za-z0-9\'\"]+)\s+\((....)\)\s*$/) {
                    $d1large_coef->[$id1] = $2;
                    $d1large_i->[$id1] = "$3$4";
                    $d1large_j->[$id1] = "$5$6";
                    $d1large_a->[$id1] = "$7$8";
                    $d1large_b->[$id1] = "$9$10";
                    $d1large_spin->[$id1] = $11;
                    $id1 = $id1 + 1;
                }
                else {
                    last;
                }
            }
            $self->{"d1large_coef"} = $d1large_coef;
            $self->{"d1large_i"} = $d1large_i;
            $self->{"d1large_j"} = $d1large_j;
            $self->{"d1large_a"} = $d1large_a;
            $self->{"d1large_b"} = $d1large_b;
            $self->{"d1large_spin"} = $d1large_spin;
        }
        elsif ($wante && /^\s*Natural\s+Population\s+Analysis:\s*$/) {
            my $npacharge = [];
            my $npashellpop = [];
            <$out>;
            my $iatom = 0;
            while (<$out>) {
                if (/^\s*\d+\s+[A-Za-z]+\s+$fltrx\s+$fltrx\s*$/) {
                    $npacharge->[$iatom] = $1;
                    $npashellpop->[$iatom] = [ $2 ];
                }
                elsif (/^\s*\d+\s+[A-Za-z]+\s+$fltrx\s+$fltrx\s+$fltrx\s*$/) {
                    $npacharge->[$iatom] = $1;
                    $npashellpop->[$iatom] = [ $2, $3 ];
                }
                elsif (/^\s*\d+\s+[A-Za-z]+\s+$fltrx\s+$fltrx\s+$fltrx\s+$fltrx\s*$/) {
                    $npacharge->[$iatom] = $1;
                    $npashellpop->[$iatom] = [ $2, $3, $4 ];
                }
                elsif (/^\s*\d+\s+[A-Za-z]+\s+$fltrx\s+$fltrx\s+$fltrx\s+$fltrx\s+$fltrx\s*$/) {
                    $npacharge->[$iatom] = $1;
                    $npashellpop->[$iatom] = [ $2, $3, $4, $5 ];
                }
                elsif (/^\s*\d+\s+[A-Za-z]+\s+$fltrx\s+$fltrx\s+$fltrx\s+$fltrx\s+$fltrx\s+$fltrx\s*$/) {
                    $npacharge->[$iatom] = $1;
                    $npashellpop->[$iatom] = [ $2, $3, $4, $5, $6 ];
                }
                elsif (/^\s*\d+\s+[A-Za-z]+\s+$fltrx\s+$fltrx\s+$fltrx\s+$fltrx\s+$fltrx\s+$fltrx\s+$fltrx\s*$/) {
                    $npacharge->[$iatom] = $1;
                    $npashellpop->[$iatom] = [ $2, $3, $4, $5, $6, $7 ];
                }
                elsif (/^\s*$/) {
                    last;
                }
                else {
                    die "AM too high to read NPA shell populations (line=$_)";
                }
                $iatom = $iatom + 1;
            }
            $self->{"npacharge"} = $npacharge;
            $self->{"npashellpop"} = $npashellpop;
        }
        elsif (/The optimization has converged/) {
            $optconverged = 1;
        }
        elsif (/^\s*Beginning displacement/) {
            # don't grab energies for displaced geometries
            $wante = 0;
        }
        elsif ($optconverged
               && /^\s+n\s+atom\s+label\s+x\s+y\s+z\s+mass\s*$/) {
            # old style geometry
            my $molstr = "";
            while (<$out>) {
                s/^ *[0-9]+:// if ($have_nodenum);
                if (! /^\s+\d+\s+(\w+)\s+$fltrx\s+$fltrx\s+$fltrx\s+/
                 && ! /^\s+\d+\s+(\w+)\s+\S+\s+$fltrx\s+$fltrx\s+$fltrx\s+/) {
                    last;
                }
                $molstr = sprintf "%s %s %16.14f %16.14f %16.14f\n",
                    ${molstr}, $1, $2 * $to_angstrom,
                    $3 * $to_angstrom, $4 * $to_angstrom;
            }
            $molecule = new Molecule($molstr);
        }
        elsif (/^\s*molecule<Molecule>:/) {
            $geometry_conversion = $to_angstrom;
        }
        elsif (/^\s*unit = \"angstrom\"\s*/) { # "
            $geometry_conversion = 1.0;
        }
        elsif ($optconverged
               && /n\s+atoms\s+(atom_labels\s+)?geometry/) {
            # new style geometry
            my $molstr = "";
            while (<$out>) {
                s/^ *[0-9]+:// if ($have_nodenum);
                if (! /^\s+\d+\s+(\w+)\s+\[\s*$fltrx\s+$fltrx\s+$fltrx\s*\]/
                 && ! /^\s+\d+\s+(\w+)\s+\"[^\"]*\"+\s+\[\s*$fltrx\s+$fltrx\s+$fltrx\s*\]/) { # " (unconfuse emacs)
                    last;
                }
                $molstr = sprintf "%s %s %16.14f %16.14f %16.14f\n",
                    ${molstr}, $1, $2 * $geometry_conversion,
                    $3 * $geometry_conversion, $4 * $geometry_conversion;
            }
            $molecule = new Molecule($molstr);
        }
        elsif (/^\s+Total (MP2 )?[Gg]radient/) {
            $state = "read grad";
            $grad = [];
            $ngrad = 0;
        }
        elsif (/^\s+Frequencies .*:\s*$/) {
            # read the frequencies
            while (<$out>) {
                s/^ *[0-9]+:// if ($have_nodenum);
                if (/^\s+\d+\s+$fltrx\s*$/) {
                    $freq->[$ifreq] = $1;
                    $ifreq++;
                }
                elsif (/THERMODYNAMIC ANALYSIS/) {
                    last;
                }
            }
            $havefreq = 1;
        }
    }
    $self->{"t1norm"} = $t1norm;
    $self->{"t1matrix2norm"} = $t12norm;
    $self->{"s2norm"} = $s2norm;
    $self->{"ccsdenergy"} = $ccsd_energy;
    $self->{"ccsd_tenergy"} = $ccsd_t_energy;
    $self->{"scfenergy"} = $scfenergy;
    $self->{"opt1energy"} = $opt1energy;
    $self->{"opt2energy"} = $opt2energy;
    $self->{"zapt2energy"} = $zapt2energy;
    $self->{"mp2energy"} = $mp2energy;
    $self->{"optconverged"} = $optconverged;
    $self->{"molecularenergy"} = $molecularenergy;
    if ($optconverged) {
        $self->{"optmolecule"} = $molecule;
    }
    $self->{"have_frequencies"} = $havefreq;
    if ($havefreq) {
        my @tmp = sort(@$freq);
        $freq = \@tmp;
        $self->{"freq"} = $freq;
    }
    my $qcinput = $self->{"qcinput"};
    my $method = $qcinput->method();
    #printf "qcinput ok = %d\n", $qcinput->ok();
    if ($method eq "MP2") {
        if ($mp2energy ne "") {
            $self->{"energy"} = $mp2energy;
        }
        elsif ($qcinput->mult() == 1) {
            $self->{"energy"} = $zapt2energy;
        }
    }
    elsif ($method eq "OPT1[2]") {
        $self->{"energy"} = $opt1energy;
    }
    elsif ($method eq "OPT2[2]") {
        $self->{"energy"} = $opt2energy;
    }
    elsif ($method eq "CCSD") {
        $self->{"energy"} = $ccsd_energy;
    }
    elsif ($method eq "CCSD(T)") {
        $self->{"energy"} = $ccsd_t_energy;
    }
    elsif ($qcinput->ok()
           && ($method eq "SCF" || $method eq "ROSCF")) {
        $self->{"energy"} = $scfenergy;
    }

    if ($self->{"energy"} eq "") {
        $self->{"energy"} = $self->{"molecularenergy"};
    }

    $self->{"ok"} = 0;
    if ($self->{"energy"} ne "") {
        if ($qcinput->optimize()) {
            if ($self->{"optconverged"}) {
                $self->{"ok"} = 1;
            }
            else {
                $error = "$error geom not converged\n";
            }
        }
        else {
            $self->{"ok"} = 1;
        }
        if ($qcinput->frequencies() && ! $havefreq) {
            $self->{"ok"} = 0;
            $error = "$error no freq\n";
        }
    }
    else {
        $error = "$error no energy\n";
    }

    if (! $qcinput->ok()) {
        $self->{"ok"} = 0;
        $error = "$error qcinput error: $qcinput->{'error'}";
    }

    $self->{"error"} = $error;
}

sub ok {
    my $self = shift;
    $self->{"ok"}
}

sub error {
    my $self = shift;
    $self->{"error"};
}

sub exists {
    my $self = shift;
    $self->{"exists"}
}

sub input {
    my $self = shift;
    $self->{"qcinput"}
}

sub inputok {
    my $self = shift;
    $self->{"qcinput"}->ok();
}

sub energy {
    my $self = shift;
    $self->{"energy"}
}

sub s2norm {
    my $self = shift;
    $self->{"s2norm"}
}

sub s2matrix1norm {
    my $self = shift;
    $self->{"s2matrix1norm"}
}

sub degenerate {
    my $self = shift;
    $self->{"degenerate"}
}

sub d1mp2 {
    my $self = shift;
    $self->{"d1mp2"}
}

sub d2mp1 {
    my $self = shift;
    $self->{"d2mp1"}
}

sub s2matrixinfnorm {
    my $self = shift;
    $self->{"s2matrixinfnorm"}
}

sub t1norm {
    my $self = shift;
    $self->{"t1norm"}
}

sub t1matrix2norm {
    my $self = shift;
    $self->{"t1matrix2norm"}
}

sub optmolecule {
    my $self = shift;
    $self->{"optmolecule"};
}

sub gradient {
    my $self = shift;
    $self->{"grad"};
}

sub frequencies {
    my $self = shift;
    $self->{"freq"}
}

sub d1large_coef {
    my $self = shift;
    $self->{"d1large_coef"}
}

sub d1large_i {
    my $self = shift;
    $self->{"d1large_i"}
}

sub d1large_j {
    my $self = shift;
    $self->{"d1large_j"}
}

sub d1large_a {
    my $self = shift;
    $self->{"d1large_a"}
}

sub d1large_b {
    my $self = shift;
    $self->{"d1large_b"}
}

sub d1large_spin {
    my $self = shift;
    $self->{"d1large_spin"}
}

sub s2large_coef {
    my $self = shift;
    $self->{"s2large_coef"}
}

sub s2large_i {
    my $self = shift;
    $self->{"s2large_i"}
}

sub s2large_a {
    my $self = shift;
    $self->{"s2large_a"}
}

sub npacharge {
    my $self = shift;
    $self->{"npacharge"}
}

sub npashellpop {
    my $self = shift;
    $self->{"npashellpop"}
}

1;
