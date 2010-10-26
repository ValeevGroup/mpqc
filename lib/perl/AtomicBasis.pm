#
eval 'exec perl $0 $*'
    if 0;

package AtomicBasis;

$fltrx = "((?:-?\\d+|-?\\d+\\.\\d*|-?\\d*\\.\\d+)(?:[eEdD][+-]?\\d+)?)";

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

    $self->{"pure_d"} = 1;
    $self->{"pure_f_plus"} = 1;
    $self->{"exponents"} = [];
    $self->{"coefficients"} = [];
    $self->{"am"} = [];
    $self->{"pure"} = [];
}

sub default_pure {
    my $self = shift;
    my $am = shift;
    if ($am == 2 && $self->{"pure_d"}) { return 1; }
    elsif ($am > 2 && $self->{"pure_f_plus"}) { return 1; }
    return 0;
}

my $amtypes = "SPDFGHIKLMN";

sub am_string_to_number {
    my $amstr = uc(shift);
    if (length($amstr) != 1) {
        die "invalid am string \"$amstr\"";
    }
    my $index = index($amtypes,$amstr);
    if ($index == -1) {
        die "am string \"$amstr\" not found";
    }
    return $index;
}

sub am_number_to_string {
    my $am = shift;
    return substr($amtypes,$am,1);
}

sub stdflt {
    my $num = shift;
    $num =~ s/D/e/;
    $num =~ s/d/e/;
    $num =~ s/E/e/;
    return $num;
}

sub read_gaussian {
    my $self = shift;
    my $file = shift;
    my $ishell = 0;
    while (<$file>) {
        if (/\*\*\*\*/) {
            last;
        }
        elsif (/([A-Za-z]+) +([0-9]+) +([0-9.]+)/) {
            my $amstr = uc($1);
            my $nprim = $2;
            my $scale = $3;
            if ($scale != 1) {
                die "cannot handle scale = $scale (must be 1)";
            }
            if ($amstr eq "SP") {
                $self->{"am"}->[$ishell]->[0] = 0;
                $self->{"am"}->[$ishell]->[1] = 1;
                $self->{"pure"}->[$ishell]->[0] = 0;
                $self->{"pure"}->[$ishell]->[1] = 0;
                foreach my $i (0..$nprim-1) {
                    while (<$file>) { last; }
                    if (/$fltrx\s+$fltrx\s+$fltrx\s*$/) {
                        my $exp = $1;
                        my $coefs = $2;
                        my $coefp = $3;
                        $self->{"exponents"}->[$ishell]->[$i]
                            = stdflt($exp);
                        $self->{"coefficients"}->[$ishell]->[$i]->[0]
                            = stdflt($coefs);
                        $self->{"coefficients"}->[$ishell]->[$i]->[1]
                            = stdflt($coefp);
                    }
                    else {
                        die "bad exponent coefficient line";
                    }
                }
            }
            else {
                my $am = am_string_to_number($amstr);
                $self->{"am"}->[$ishell]->[0] = $am;
                $self->{"pure"}->[$ishell]->[0] = $self->default_pure($am);
                foreach my $i (0..$nprim-1) {
                    while (<$file>) { last; }
                    if (/$fltrx\s+$fltrx\s*$/) {
                        my $exp = $1;
                        my $coef = $2;
                        $self->{"exponents"}->[$ishell]->[$i]
                            = stdflt($exp);
                        $self->{"coefficients"}->[$ishell]->[$i]->[0]
                            = stdflt($coef);
                    }
                    else {
                        die "bad exponent coefficient line";
                    }
                }
            }
            $ishell++;
        }
        else {
            die "could not parse line $_";
        }
    }
}

sub nshell {
    my $self = shift;
    my @exp = @{$self->{"exponents"}};
    return $#exp + 1;
}

sub nprim {
    my $self = shift;
    my $ishell = shift;
    my @exp = @{$self->{"exponents"}->[$ishell]};
    return $#exp + 1;
}

sub ncon {
    my $self = shift;
    my $ishell = shift;
    my @exp = @{$self->{"am"}->[$ishell]};
    return $#exp + 1;
}

sub exp {
    my $self = shift;
    my $ishell = shift;
    my $iprim = shift;
    my $exp = $self->{"exponents"}->[$ishell]->[$iprim];
    return $exp;
}

sub coef {
    my $self = shift;
    my $ishell = shift;
    my $iprim = shift;
    my $icon = shift;
    my $coef = $self->{"coefficients"}->[$ishell]->[$iprim]->[$icon];
    return $coef;
}

sub am {
    my $self = shift;
    my $ishell = shift;
    my $icon = shift;
    my $am = $self->{"am"}->[$ishell]->[$icon];
    return $am;
}

sub pure {
    my $self = shift;
    my $ishell = shift;
    my $icon = shift;
    my $pure = $self->{"pure"}->[$ishell]->[$icon];
    return $pure;
}

sub amstr {
    my $self = shift;
    my $ishell = shift;
    my $icon = shift;
    my $am = $self->{"am"}->[$ishell]->[$icon];
    return am_number_to_string($am);
}

sub write_keyval {
    my $self = shift;
    my $file = shift;
    foreach my $ishell (0..$self->nshell()-1) {
        print $file "  (type: [";
        # write out am (and puream)
        foreach my $icon (0..$self->ncon($ishell)-1) {
            if ($icon > 0) { print $file " "; }
            if ($self->pure($ishell,$icon)) {
                printf $file "(am=%s puream=1)",
                lc($self->amstr($ishell,$icon));
            }
            else {
                printf $file "am=%s", lc($self->amstr($ishell,$icon));
            }
        }
        print $file "]\n";
        print $file "   {exp";
        # write out coef:0...
        foreach my $icon (0..$self->ncon($ishell)-1) {
            print $file " coef:$icon";
        }
        print "} = {\n";
        foreach my $iprim (0..$self->nprim($ishell)-1) {
            printf $file "   %s", $self->exp($ishell,$iprim);
            # write out coefficients
            foreach my $icon (0..$self->ncon($ishell)-1) {
                printf $file " %s", $self->coef($ishell,$iprim,$icon);
            }
            printf $file "\n";
        }
        print $file "  })\n";
    }
}

1;
