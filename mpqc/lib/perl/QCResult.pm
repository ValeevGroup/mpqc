#
eval 'exec perl $0 $*'
    if 0;

require QCParse;
require Molecule;

##########################################################################

package QCResult;

# this seems to not work as expected
$fltrx = "([-+]?(?:\\d+\\.\\d*|\\d+|\\.\\d+)(?:[eEdD][+-]?\\d+)?)";

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
        while (<OUTFILE>) {
            if (/^\s*MPQC:/) {
                $self->parse_mpqc(\*OUTFILE);
                last;
            }
            elsif (/^\s*Entering Gaussian System/) {
                $self->parse_g94(\*OUTFILE);
                last;
            }
        }
        close(<OUTFILE>);
    }
}

sub parse_g94 {
    my $self = shift;
    my $out = shift;
    my $scfenergy = "";
    my $mp2energy = "";
    my $optconverged = 0;
    my $molecule;
    my $havefreq = 0;
    my $freq = [];
    my $ifreq = 0;
    while (<$out>) {
        if (/^\s*SCF Done:  E\(RHF\) =\s*$fltrx\s/) {
            $scfenergy = $1;
        }
        elsif (/^\s*E2\s*=\s*$fltrx\s*EUMP2\s*=\s*$fltrx/) {
            $mp2energy = $2;
        }
        elsif (/^\s*-- Stationary point found./) {
            $optconverged = 1;
        }
        elsif ($optconverged && /^\s*Input orientation:/) {
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
    elsif ($method eq "SCF"
           || $method eq "ROSCF") {
        $self->{"energy"} = $scfenergy;
    }

    $self->{"ok"} = 0;
    if ($self->{"energy"} ne "") {
        if ($qcinput->optimize()) {
            $self->{"ok"} = 1 if $self->{"optconverged"};
        }
        else {
            $self->{"ok"} = 1;
        }
        if ($qcinput->frequencies() && ! $havefreq) {
            $self->{"ok"} = 0;
        }
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
    while (<$out>) {
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
        elsif ($wante && /Value of the MolecularEnergy:\s+$fltrx/) {
            $molecularenergy = $1;
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
            my $molstr = "";
            while (<$out>) {
                if (! /^\s+\d+\s+(\w+)\s+$fltrx\s+$fltrx\s+$fltrx\s+/
                 && ! /^\s+\d+\s+(\w+)\s+\S+\s+$fltrx\s+$fltrx\s+$fltrx\s+/) {
                    last;
                }
                $molstr = "${molstr} $1 $2 $3 $4\n";
            }
            $molecule = new Molecule($molstr);
        }
        elsif (/^\s+Total (MP2 )?[Gg]radient/
               || /^\s*Gradient of the MolecularEnergy:/) {
            $state = "read grad";
            $grad = [];
            $ngrad = 0;
        }
        elsif (/^\s+Frequencies .*:\s*$/) {
            # read the irrep
            $_ = <$out>;
            if (! /\s*A\s*/) {
                print "More than one freq irrep: need to modify QCResult.pm\n";
                exit 1;
            }
            $freq = [];
            # read the frequencies
            while (<$out>) {
                if (/^\s+\d+\s+$fltrx\s*$/) {
                    $freq->[$ifreq] = $1;
                    $ifreq++;
                }
                else {
                    last;
                }
            }
            $havefreq = 1;
        }
    }
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
            $self->{"ok"} = 1 if $self->{"optconverged"};
        }
        else {
            $self->{"ok"} = 1;
        }
        if ($qcinput->frequencies() && ! $havefreq) {
            $self->{"ok"} = 0;
        }
    }
}

sub ok {
    my $self = shift;
    $self->{"ok"}
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
