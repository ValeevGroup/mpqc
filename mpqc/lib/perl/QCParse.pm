#
eval 'exec perl $0 $*'
    if 0;

require Molecule;

##########################################################################

package QCParse;
$debug = 0;

sub testparse {
    my $parse = new QCParse;

    my $string = "x:
 xval
test_basis: STO-3G 6-311G**
charge: 1
method: scf
basis: sto-3g
state: 3b1
molecule:
  H 0 0.0000001 1.00000001
  H 0 0 -1
gradient: yes
optimize: no
frequencies: yes
properties: NPA
y:
yval
z: zval1 zval2
zval3
h:
0 a
1
 2  c";

    print "string:\n--------------\n$string\n--------------\n";

    $parse->parse_string($string);
    $parse->doprint();

    my @t = $parse->value_as_array('h');
    print "-----------------\n";
    for ($i = 0; $i <= $#t; $i++) {
        print "$i: $t[$i]\n";
    }
    print "-----------------\n";

    @t = $parse->value_as_lines('h');
    print "-----------------\n";
    for ($i = 0; $i <= $#t; $i++) {
        print "$i: $t[$i]\n";
    }
    print "-----------------\n";

    my $qcinp = new QCInput($parse);
    my $test_basis = $parse->value("test_basis");
    my @test_basis_a = $parse->value_as_array("test_basis");
    my $state = $qcinp->state();
    my $mult = $qcinp->mult();
    my $method = $qcinp->method();
    my $charge = $qcinp->charge();
    my $basis = $qcinp->basis();
    my $gradient = $qcinp->gradient();
    my $frequencies = $qcinp->frequencies();
    my $optimize = $qcinp->optimize();
    my $natom = $qcinp->n_atom();
    foreach $i (@test_basis_a) {
        print "test_basis_a: $i\n";
    }
    print "test_basis = $test_basis\n";
    print "state = $state\n";
    print "mult = $mult\n";
    print "method = $method\n";
    print "basis = $basis\n";
    print "optimize = $optimize\n";
    print "gradient = $gradient\n";
    print "frequencies = $frequencies\n";
    print "natom = $natom\n";
    for ($i = 0; $i < $natom; $i++) {
        printf "%s %14.8f %14.8f %14.8f\n", $qcinp->element($i),
                                $qcinp->position($i,0),
                                $qcinp->position($i,1),
                                $qcinp->position($i,2);
    }
    printf "qcinp errors: %s\n", $qcinp->error();

    my $inpwr = new MPQCInputWriter($qcinp);
    printf "MPQC input:\n%s", $inpwr->input_string();
}

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {};
    bless $self, $class;
    $self->initialize();
    return $self;
}

sub initialize {
    my $self = shift;
    $self->{'keyval'} = {};
    $self->{'error'} = "";
}

sub parse_file {
    my $self = shift;
    my $file = shift;
    if (! -f "$file") {
        $self->{"ok"} = 0;
        $self->error("File $file not found.");
        return;
    }
    open(INPUT, "<$file");
    my $string = "";
    while (<INPUT>) {
        $string = "$string$_";
    }
    close(INPUT);
    #print "Got file:\n$string\n";
    $self->parse_string($string);
    $self->{"ok"} = 1;
}

sub write_file {
    my $self = shift;
    my $file = shift;
    my $keyval = $self->{'keyval'};
    my @keys = keys(%$keyval);
    open(OUTPUT, ">$file");
    foreach $key (@keys) {
        my $value = $keyval->{$key};
        print OUTPUT "${key}:\n";
        print OUTPUT "$value\n";
    }
    close(OUTPUT);
}

sub parse_string {
    my $self = shift;
    my $string = shift;
    my $value = "";
    my $keyword = "";
    $string = "$string\n";
    while ($string) {
        $string =~ s/^[^\n]*\n//;
        $_ = $&;
        s/#.*//;
        if (/^\s*(\w+)\s*:\s*(.*)\s*$/) {
            $self->add($keyword, $value);
            $keyword = $1;
            $value = $2;
        }
        elsif (/^\s*$/) {
            $self->add($keyword, $value);
            $keyword = "";
            $value = "";
        }
        else {
            $value = "$value$_";
        }
    }
    $self->add($keyword, $value);
}

sub add {
    my $self = shift;
    my $keyword = shift;
    my $value = shift;
    if ($keyword ne "") {
        $self->{'keyval'}{$keyword} = $value;
        printf("%s = %s\n", $keyword, $value) if ($debug);
    }
}

# returns the value of the keyword
sub value {
    my $self = shift;
    my $keyword = shift;
    my $keyval = $self->{'keyval'};
    my $value = $keyval->{$keyword};
    return $value;
}

# sets the value of the keyword
sub set_value {
    my $self = shift;
    my $keyword = shift;
    my $value = shift;
    my $keyval = $self->{'keyval'};
    $keyval->{$keyword} = $value;
    return $value;
}

# returns the value of the keyword
sub boolean_value {
    my $self = shift;
    my $keyword = shift;
    my $keyval = $self->{'keyval'};
    $_ = $keyval->{$keyword};
    return "1" if (/^\s*(y|yes|1|true|t)\s*$/i);
    return "0" if (/^\s*(n|no|0|false|f|)\s*$/i);
    "";
}

# returns an array of whitespace delimited tokens
sub value_as_array {
    my $self = shift;
    my $keyword = shift;
    my $keyval = $self->{'keyval'};
    my $value = $keyval->{$keyword};
    my @array = ();
    $i = 0;
    $value =~ s/^\s+$//;
    while ($value ne '') {
        $value =~ s/^\s*(\S+)\s*//s;
        $array[$i] = $1;
        $i++;
    }
    return @array;
}

# returns an array reference of whitespace delimited tokens
sub value_as_arrayref {
    my $self = shift;
    my $keyword = shift;
    my $keyval = $self->{'keyval'};
    my $value = $keyval->{$keyword};
    my $array = [];
    $i = 0;
    $value =~ s/^\s+$//;
    while ($value ne '') {
        $value =~ s/^\s*(\S+)\s*//s;
        $array->[$i] = $1;
        $i++;
    }
    return $array;
}

# returns an array of lines
sub value_as_lines {
    my $self = shift;
    my $keyword = shift;
    my $keyval = $self->{'keyval'};
    my $value = $keyval->{$keyword};
    my @array = ();
    $i = 0;
    while ($value) {
        $value =~ s/^\s*(.*)\s*\n//;
        $array[$i] = $1;
        $i++;
    }
    return @array;
}

# returns 1 if the input file existed
sub ok {
    my $self = shift;
    $self->{"ok"};
}

sub display {
    my $self = shift;
    my @keys = @_ ? @_ : sort keys %$self;
    foreach $key (@keys) {
        print "\t$key => $self->{$key}\n";
    }
}

sub doprint {
    my $self = shift;
    print "QCParse:\n";
    my $keyval = $self->{'keyval'};
    foreach $i (keys %$keyval) {
        my $val = $keyval->{$i};
        $val =~ s/\n/\\n/g;
        print "keyword = $i, value = $val\n";
    }
}

sub error {
    my $self = shift;
    my $msg = shift;
    $self->{"error"} = "$self->{'error'}$msg";
}

##########################################################################

package QCInput;
$debug = 0;

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
    my $parser = shift;
    if ($parser eq "") {
        $parser = new QCParse;
    }
    $self->{"parser"} = $parser;
    $self->{"error"} = $parser->error();

    $self->{"molecule"} = new Molecule($parser->value("molecule"));
}

sub error {
    my $self = shift;
    my $msg = shift;
    $self->{"error"} = "$self->{'error'}$msg";
}

sub checkpoint {
    my $self = shift;
    my $bval = $self->{"parser"}->boolean_value("checkpoint");
    my $val = $self->{"parser"}->value("checkpoint");
    if ($val ne "" && $bval eq "") {
        $self->error("Bad value for checkpoint: $val");
        $bval = "0";
    }
    elsif ($val eq "") {
        $bval = "1";
    }
    $bval;
}

sub restart {
    my $self = shift;
    my $bval = $self->{"parser"}->boolean_value("restart");
    my $val = $self->{"parser"}->value("restart");
    if ($val ne "" && $bval eq "") {
        $self->error("Bad value for restart: $val");
        $bval = "0";
    }
    elsif ($val eq "") {
        $bval = "1";
    }
    $bval;
}

sub label {
    my $self = shift;
    $self->{"parser"}->value("label");
}

sub charge {
    my $self = shift;
    $_ = $self->{"parser"}->value("charge");
    s/^\s+//;
    s/\s+$//;
    s/^\+//;
    if (/^$/) { $_ = "0"; }
    if (! /^-?\d+$/) {
        $self->error("Bad charge: $_ (using 0)\n");
        $_ = "0";
    }
    $_;
}

sub method {
    my $self = shift;
    $_ = $self->{"parser"}->value("method");
    s/^\s+//;
    s/\s+$//;
    if ($_ eq "") {
        $self->error("No method given (using default).\n");
        $_ = "SCF";
    }
    tr/a-z/A-Z/;
    $_;
}

sub symmetry {
    my $self = shift;
    $_ = $self->{"parser"}->value("symmetry");
    s/^\s*//;
    s/\s*$//;
    uc $_;
}

sub memory {
    my $self = shift;
    $_ = $self->{"parser"}->value("memory");
    s/^\s*//;
    s/\s*$//;
    if ($_ eq "") {
        $_ = 32000000;
    }
    $_;
}

sub state {
    my $self = shift;
    $_ = $self->{"parser"}->value("state");
    s/^\s*//;
    s/\s*$//;
    uc $_;
}

sub mult {
    my $self = shift;
    $_ = $self->state();
    s/^\s*(\d+)/\1/;
    if (/^\s*$/) {
        $_ = 1;
    }
    $_;
}

sub basis {
    my $self = shift;
    $_ = $self->{"parser"}->value("basis");
    s/^\s+//;
    s/\s+$//;
    if ($_ eq "") {
        $self->error("No basis given (using default).\n");
        $_ = "STO-3G";
    }
    $_;
}

sub auxbasis {
    my $self = shift;
    $_ = $self->{"parser"}->value("auxbasis");
    s/^\s+//;
    s/\s+$//;
    if ($_ eq "") {
        $self->error("No auxiliary basis given (using default).\n");
        $_ = "STO-3G";
    }
    $_;
}

sub grid {
    my $self = shift;
    $_ = $self->{"parser"}->value("grid");
    s/^\s+//;
    s/\s+$//;
    if ($_ eq "") {
        $_ = "default";
    }
    $_;
}

sub gradient {
    my $self = shift;
    my $bval = $self->{"parser"}->boolean_value("gradient");
    if ($bval eq "") {
        my $val = $self->{"parser"}->value("gradient");
        $self->error("Bad value for gradient: $val");
    }
    $bval;
}

sub fzc {
    my $self = shift;
    $_ = $self->{"parser"}->value("fzc");
    s/^\s+//;
    s/\s+$//;
    if ($_ eq "") {
        $_ = 0;
    }
    $_;
}

sub fzv {
    my $self = shift;
    $_ = $self->{"parser"}->value("fzv");
    s/^\s+//;
    s/\s+$//;
    if ($_ eq "") {
        $_ = 0;
    }
    $_;
}

sub docc {
    my $self = shift;
    $_ = $self->{"parser"}->value("docc");
    s/^\s+//;
    s/\s+$//;
    if ($_ eq "" || $_ eq "-") {
        $_ = "auto";
    }
    $_;
}

sub socc {
    my $self = shift;
    $_ = $self->{"parser"}->value("socc");
    s/^\s+//;
    s/\s+$//;
    if ($_ eq "" || $_ eq "-") {
        $_ = "auto";
    }
    $_;
}

sub optimize {
    my $self = shift;
    my $bval = $self->{"parser"}->boolean_value("optimize");
    if ($bval eq "") {
        my $val = $self->{"parser"}->value("optimize");
        $self->error("Bad value for optimize: $val");
    }
    $bval;
}

# returns "" if orthog_method not set
sub orthog_method {
    my $self = shift;
    my $bval = $self->{"parser"}->value("orthog_method");
    $bval;
}

# returns "" if lindep_tol not set
sub lindep_tol {
    my $self = shift;
    my $bval = $self->{"parser"}->value("lindep_tol");
    $bval;
}

sub transition_state {
    my $self = shift;
    my $bval = $self->{"parser"}->boolean_value("transition_state");
    if ($bval eq "") {
        my $val = $self->{"parser"}->value("transition_state");
        $self->error("Bad value for transtion_state: $val");
    }
    $bval;
}

sub frequencies {
    my $self = shift;
    my $bval = $self->{"parser"}->boolean_value("frequencies");
    if ($bval eq "") {
        my $val = $self->{"parser"}->value("frequencies");
        $self->error("Bad value for frequencies: $val");
    }
    $bval;
}

sub axyz_lines {
    my $self = shift;
    $self->molecule()->string();
}

sub molecule() {
    my $self = shift;
    return $self->{"molecule"};
}

sub n_atom {
    my $self = shift;
    printf "QCInput: returning natom = %d\n", $self->{"natom"} if ($debug);
    $self->molecule()->n_atom();
}

sub element {
    my $self = shift;
    $self->molecule()->element(@_);
}

sub position {
    my $self = shift;
    $self->molecule()->position(@_);
}

sub write_file {
    my $self = shift;
    my $file = shift;
    my $parser = $self->{'parser'};
    $parser->write_file($file);
}

sub mode_following() {
    my $self = shift;
    return scalar($self->{"parser"}->value_as_array("followed")) != 0;
}

# returns 1 if the input file existed
sub ok {
    my $self = shift;
    $self->{"parser"}->{"ok"};
}

##########################################################################

package InputWriter;

# Input Writer is abstract
sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {};
    bless $self, $class;
    $self->initialize(@_);
    return $self;
}

sub initialize() {
    my $self = shift;
    my $qcinput = shift;
    $self->{"qcinput"} = $qcinput;
}

# this should be overridden
sub input_string() {
    "";
}

sub write_input() {
    my $self = shift;
    my $file = shift;
    my $input = $self->input_string();
    open(OUTPUT,">$file");
    printf OUTPUT "%s", $input;
    close(OUTPUT);
}

sub write_qcinput {
    my $self = shift;
    my $file = shift;
    my $qcinput = $self->{'qcinput'};
    $qcinput->write_file($file);
}

##########################################################################

package MPQCInputWriter;
@ISA = qw( InputWriter );
%methodmap = ("MP2-R12/A" => "MBPT2_R12",
              "MP2-R12/A'" => "MBPT2_R12",
              "MP2" => "MBPT2",
              "OPT1[2]" => "MBPT2",
              "OPT2[2]" => "MBPT2",
              "ZAPT2" => "MBPT2",
              "MP2V1" => "MBPT2",
              "OPT1[2]V1" => "MBPT2",
              "OPT2[2]V1" => "MBPT2",
              "ZAPT2V1" => "MBPT2",
              "MP2V2" => "MBPT2",
              "OPT1[2]V2" => "MBPT2",
              "OPT2[2]V2" => "MBPT2",
              "ZAPT2V2" => "MBPT2",
              "MP2V2LB" => "MBPT2",
              "OPT1[2]V2LB" => "MBPT2",
              "OPT2[2]V2LB" => "MBPT2",
              "ZAPT2V2LB" => "MBPT2",
              "ROSCF" => "SCF",
              "SCF" => "SCF",
              "UHF" => "UHF",
              "CLHF" => "CLHF",
              "HSOSHF" => "HSOSHF",
              "HF" => "SCF",
              "HFK" => "DFT",
              "XALPHA" => "DFT",
              "HFS"    => "DFT",
              "HFB"    => "DFT",
              "HFG96"  => "DFT",
              "BLYP"   => "DFT",
              "B3LYP"  => "DFT",
              "B3PW91"  => "DFT",
              "PBE"    => "DFT",
              "PW91"    => "DFT",
              "SPZ81"   => "DFT",
              "B3P86"   => "DFT",
              "BP86"    => "DFT",
              "BPW91"   => "DFT",
              "CLHFK" => "DFT",
              "CLXALPHA" => "DFT",
              "CLHFS"    => "DFT",
              "CLHFB"    => "DFT",
              "CLHFG96"  => "DFT",
              "CLBLYP"   => "DFT",
              "CLB3LYP"  => "DFT",
              "CLB3PW91"  => "DFT",
              "CLPBE"    => "DFT",
              "CLPW91"    => "DFT",
              "SPZ81"   => "DFT",
              "B3P86"   => "DFT",
              "BP86"    => "DFT",
              "BPW91"   => "DFT",
              "HSOSHFK" => "DFT",
              "HSOSXALPHA" => "DFT",
              "HSOSHFS"    => "DFT",
              "HSOSHFB"    => "DFT",
              "HSOSHFG96"  => "DFT",
              "HSOSBLYP"   => "DFT",
              "HSOSB3LYP"  => "DFT",
              "HSOSB3PW91"  => "DFT",
              "HSOSPBE"    => "DFT",
              "HSOSPW91"    => "DFT",
              "HSOSSPZ81"   => "DFT",
              "HSOSB3P86"   => "DFT",
              "HSOSBP86"    => "DFT",
              "HSOSBPW91"   => "DFT",
              "UHFK" => "DFT",
              "UXALPHA" => "DFT",
              "UHFS"    => "DFT",
              "UHFB"    => "DFT",
              "UHFG96"  => "DFT",
              "UBLYP"   => "DFT",
              "UB3LYP"  => "DFT",
              "UB3PW91"  => "DFT",
              "UPBE"    => "DFT",
              "UPW91"    => "DFT",
              "USPZ81"   => "DFT",
              "UB3P86"   => "DFT",
              "UBP86"    => "DFT",
              "UBPW91"   => "DFT",
             );
%mbpt2r12stdapproxmap = ("MP2-R12/A" => "A",
                         "MP2-R12/A'" => "A'",
                         );
%mbpt2map = ("MP2" => "mp",
             "OPT1[2]" => "opt1",
             "OPT2[2]" => "opt2",
             "ZAPT2" => "zapt",
             "MP2V1" => "mp",
             "OPT1[2]V1" => "opt1",
             "OPT2[2]V1" => "opt2",
             "ZAPT2V1" => "zapt",
             "MP2V2" => "mp",
             "OPT1[2]V2" => "opt1",
             "OPT2[2]V2" => "opt2",
             "ZAPT2V2" => "zapt",
             "MP2V2LB" => "mp",
             "OPT1[2]V2LB" => "opt1",
             "OPT2[2]V2LB" => "opt2",
             "ZAPT2V2LB" => "zapt");
%mbpt2algmap = ("MP2" => "",
                "OPT1[2]" => "",
                "OPT2[2]" => "",
                "ZAPT2" => "",
                "MP2V1" => "v1",
                "OPT1[2]V1" => "v1",
                "OPT2[2]V1" => "v1",
                "ZAPT2V1" => "v1",
                "MP2V2" => "v2",
                "OPT1[2]V2" => "v2",
                "OPT2[2]V2" => "v2",
                "ZAPT2V2" => "v2",
                "MP2V2LB" => "v2lb",
                "OPT1[2]V2LB" => "v2lb",
                "OPT2[2]V2LB" => "v2lb",
                "ZAPT2V2LB" => "v2lb");
$debug = 0;

sub new {
    my $this = shift;
    my $class = ref($this) || $this;
    my $self = {};
    bless $self, $class;

    $self->initialize(@_);
    return $self;
}

sub initialize() {
    my $self = shift;
    my $qcinput = shift;
    $self->{"qcinput"} = $qcinput;
}

sub docc_string() {
    my $self = shift;
    my $qcinput = $self->{"qcinput"};
    my $occs = $qcinput->docc();
    if ($occs eq "auto") { return ""; }
    $occs =~ s/,/ /g;
    "docc = [ $occs ]";
}

sub socc_string() {
    my $self = shift;
    my $qcinput = $self->{"qcinput"};
    my $occs = $qcinput->socc();
    if ($occs eq "auto") { return ""; }
    $occs =~ s/,/ /g;
    "socc = [ $occs ]";
}

sub input_string() {
    my $self = shift;
    my $qcinput = $self->{"qcinput"};
    my $qcparse = $qcinput->{"parser"};

    my $use_cints = 0;
    my $do_cca = $qcparse->value("do_cca");

    printf "molecule = %s\n", $qcparse->value("molecule") if ($debug);

    my $symmetry = $qcinput->symmetry();
    my $mol = "% molecule specification";
    $mol = "$mol\nmolecule<Molecule>: (";
    $symmetry = lc $symmetry if ($symmetry eq "AUTO");
    if ($qcinput->frequencies()) {
        $mol = "$mol\n  symmetry = C1";
    }
    else {
        $mol = "$mol\n  symmetry = $symmetry";
    }
    $mol = "$mol\n  unit = angstrom";
    $mol = "$mol\n  { atoms geometry } = {";
    printf "MPQCInputWriter: natom = %d\n", $qcinput->n_atom() if ($debug);
    my $i;
    for ($i = 0; $i < $qcinput->n_atom(); $i++) {
        $mol = sprintf "%s\n    %2s     [ %18.12f %18.12f %18.12f ]",
                       $mol, $qcinput->element($i),
                       $qcinput->position($i,0),
                       $qcinput->position($i,1),
                       $qcinput->position($i,2);
    }
    $mol = "$mol\n  }";
    $mol = "$mol\n)\n";

    my $basis = "% basis set specification";
    $basis = "$basis\nbasis<GaussianBasisSet>: (";
    $basis = sprintf "%s\n  name = \"%s\"", $basis, $qcinput->basis();
    $basis = "$basis\n  molecule = \$:molecule";
    $basis = "$basis\n)\n";

    my $integrals = "";
    if($do_cca) {
      $integrals = "% using cca integrals";
      $integrals = "$integrals\nintegrals<IntegralCCA>: (";
      my $buffer_type = $qcparse->value("integral_buffer");
      if( $buffer_type ne "opaque" && $buffer_type ne "array" ) {
        $buffer_type = "opaque";
      }     
      my $int_package = $qcparse->value("integral_package");
      if( $int_package ne "intv3" && $int_package ne "cints" ) {
        $int_package = "intv3";
      }
      $integrals = "$integrals\n  integral_buffer = $buffer_type";
      $integrals = "$integrals\n  integral_package = $int_package";
      $integrals = "$integrals\n  evaluator_factory = MPQC.IntegralEvaluatorFactory";
      $integrals = "$integrals\n  molecule = \$:molecule";
      $integrals = "$integrals\n)\n";
    }

    my $fixed = $qcparse->value_as_arrayref("fixed");
    my $followed = $qcparse->value_as_arrayref("followed");
    if (scalar(@{$fixed}) != 0) {
        $fixed = $self->mpqc_fixed_coor($fixed);
    }
    else {
        $fixed = "";
    }
    if (scalar(@{$followed}) != 0) {
        $followed = $self->mpqc_followed_coor($followed);
    }
    else {
        $followed = "";
    }

    my $coor = "  % molecular coordinates for optimization";
    $coor = "$coor\n  coor<SymmMolecularCoor>: (";
    $coor = "$coor\n    molecule = \$:molecule";
    $coor = "$coor\n    generator<IntCoorGen>: (";
    $coor = "$coor\n      molecule = \$:molecule";
    $coor = "$coor\n    )";
    $coor = "$coor$followed";
    $coor = "$coor$fixed";
    $coor = "$coor\n  )\n";

    my $charge = $qcinput->charge();
    my $mult = $qcinput->mult();
    my $docc = $self->docc_string();
    my $socc = $self->socc_string();

    my $grid = $qcinput->grid();

    my $memory = $qcinput->memory();
    my $inputmethod = $methodmap{uc($qcinput->method())};
    my $method = "$inputmethod";
    $method = "SCF" if ($method eq "");
    my $openmethod = substr(uc($qcinput->method()),0,4);
    if (substr($openmethod,0,2) eq "CL") { $openmethod = "CL"; }
    if (substr($openmethod,0,1) eq "U") { $openmethod = "U"; }
    if ($method eq "SCF") {
        if ($openmethod eq "U") {
            $method = "UHF";
        }
        elsif ($openmethod eq "CL") {
            $method = "CLHF";
        }
        elsif ($openmethod eq "HSOS") {
            $method = "HSOSHF";
        }
        elsif ($qcinput->mult() == 1) {
            $method = "CLHF";
            $openmethod = "CL";
        }
        else {
            $method = "HSOSHF";
            $openmethod = "HSOS";
        }
    }
    my $functional;
    if ($method eq "DFT") {
        $functional = uc($qcinput->method());
        if ($openmethod eq "U") {
            $method = "UKS";
            $functional = substr($functional,1);
        }
        elsif ($openmethod eq "CL") {
            $method = "CLKS";
            $functional = substr($functional,2);
        }
        elsif ($openmethod eq "HSOS") {
            $method = "HSOSKS";
            $functional = substr($functional,4);
        }
        elsif ($qcinput->mult() == 1) {
            $method = "CLKS";
            $openmethod = "CL";
        }
        else {
            $method = "UKS";
            $openmethod = "U";
        }
    }
    my $orthog_method = $qcinput->orthog_method();
    my $lindep_tol = $qcinput->lindep_tol();
    my $mole = "  do_energy = yes";
    if ($qcinput->gradient()) {
        $mole = "$mole\n  do_gradient = yes";
    }
    else {
        $mole = "$mole\n  do_gradient = no";
    }
    if($do_cca) {
      $mole = "$mole\n  do_cca = yes";
    }
    $mole = "$mole\n  % method for computing the molecule's energy";
    $mole = "$mole\n  mole<$method>: (";
    $mole = "$mole\n    molecule = \$:molecule";
    $mole = "$mole\n    basis = \$:basis";
    $mole = "$mole\n    coor = \$..:coor";
    $mole = "$mole\n    memory = $memory";
    if($do_cca) {
      $mole = "$mole\n    integrals = \$:integrals";
    }
    if ($inputmethod eq "SCF" || $inputmethod eq "UHF"
        || $method eq "CLKS" || $method eq "UKS" || $method eq "HSOSKS") {
        $mole = "$mole\n    total_charge = $charge";
        $mole = "$mole\n    multiplicity = $mult";
        $mole = "$mole\n    print_npa = yes";
        if ($docc ne "") {$mole = "$mole\n    $docc";}
        if ($socc ne "") {$mole = "$mole\n    $socc";}
        if ($orthog_method ne "" ) {
            $mole = "$mole\n    orthog_method = $orthog_method";
        }
        if ($lindep_tol ne "" ) {
            $mole = "$mole\n    lindep_tol = $lindep_tol";
        }
    }
    if ($method eq "CLKS" || $method eq "UKS" || $method eq "HSOSKS") {
        $mole = "$mole\n    functional<StdDenFunctional>: name = \"$functional\"";
    }
    if (($method eq "CLKS" || $method eq "UKS" || $method eq "HSOSKS")
        && $grid ne "default") {
        $mole = "$mole\n    integrator<RadialAngularIntegrator>: (grid = $grid)";
    }
    if ($method eq "MBPT2_R12") {
        my $stdapprox = $mbpt2r12stdapproxmap{uc($qcinput->method())};
        my $auxbasis = $qcinput->auxbasis();
        my $fzc = $qcinput->fzc();

        $mole = sprintf "%s\n    stdapprox = \"%s\"", $mole, $stdapprox;
        $mole = "$mole\n    integrals<IntegralCints>: ()";
        $mole = "$mole\n    nfzc = $fzc";
        # don't write an auxbasis if the auxbasis is the same as the basis set.
        # this will speed up the calculation
        if ("$auxbasis" ne "" && "$auxbasis" ne $qcinput->basis()) {
            $mole = "$mole\n    aux_basis<GaussianBasisSet>: (";
            $mole = sprintf "%s\n      name = \"%s\"", $mole, $auxbasis;
            $mole = "$mole\n      molecule = \$:molecule";
            $mole = "$mole\n    )\n";
        }
        $mole = append_reference($mole,"CLHF",$charge,$mult,$memory,$orthog_method,
                                 $lindep_tol,$docc,$socc,"DZ (Dunning)");
        $use_cints = 1;
    }
    elsif ($method eq "MBPT2") {
        my $fzc = $qcinput->fzc();
        my $fzv = $qcinput->fzv();
        my $mbpt2method = $mbpt2map{uc($qcinput->method())};
        my $mbpt2algorithm = $mbpt2algmap{uc($qcinput->method())};
        $mole = "$mole\n    method = $mbpt2method";
        if ($mbpt2algorithm ne "") {
            $mole = "$mole\n    algorithm = $mbpt2algorithm";
        }
        $mole = "$mole\n    nfzc = $fzc";
        $mole = "$mole\n    nfzv = $fzv";
        my $refmethod = "";
        if ($qcinput->mult() == 1) {
            $refmethod = "CLHF";
        }
        else {
            $refmethod = "HSOSHF";
        }
        $mole = append_reference($mole,$refmethod,$charge,$mult,$memory,$orthog_method,
                                 $lindep_tol,$docc,$socc,"STO-3G");
    }
    elsif (! ($basis =~ /^STO/
              || $basis =~ /^MI/
              || $basis =~ /^\d-\d1G$/)) {
        my $guessmethod = "${openmethod}HF";
        $mole = "$mole\n    guess_wavefunction<$guessmethod>: (";
        $mole = "$mole\n      molecule = \$:molecule";
        $mole = "$mole\n      total_charge = $charge";
        $mole = "$mole\n      multiplicity = $mult";
        if ($docc ne "") {$mole = "$mole\n      $docc";}
        if ($socc ne "") {$mole = "$mole\n      $socc";}
        $mole = "$mole\n      basis<GaussianBasisSet>: (";
        $mole = "$mole\n        molecule = \$:molecule";
        $mole = "$mole\n        name = \"STO-3G\"";
        $mole = "$mole\n      )";
        $mole = "$mole\n      memory = $memory";
        if($do_cca) {
          $mole = "$mole\n      integrals = \$:integrals";
        }
        $mole = "$mole\n    )";
    }
    if ($qcinput->frequencies()) {
        $mole = "$mole\n    hessian<FinDispMolecularHessian>: (";
        if ($symmetry ne "C1") {
            $mole="$mole\n      point_group<PointGroup>: symmetry = $symmetry";
        }
        $mole = "$mole\n      checkpoint = no";
        $mole = "$mole\n      restart = no";
        $mole = "$mole\n    )";
    }
    $mole = "$mole\n  )\n";

    my $opt;
    if ($qcinput->optimize()) {
        $opt = "  optimize = yes";
    }
    else {
        $opt = "  optimize = no";
    }
    my $optclass, $updateclass;
    if ($qcinput->transition_state()) {
        $optclass = "EFCOpt";
        $updateclass = "PowellUpdate";
    }
    else {
        $optclass = "QNewtonOpt";
        $updateclass = "BFGSUpdate";
    }
    $opt = "$opt\n  % optimizer object for the molecular geometry";
    $opt = "$opt\n  opt<$optclass>: (";
    $opt = "$opt\n    max_iterations = 20";
    $opt = "$opt\n    function = \$..:mole";
    if ($qcinput->transition_state()) {
        $opt = "$opt\n    transition_state = yes";
        if ($qcinput->mode_following()) {
            $opt = "$opt\n    hessian = [ [ -0.1 ] ]";
            $opt = "$opt\n    mode_following = yes";
        }
    }
    $opt = "$opt\n    update<$updateclass>: ()";
    $opt = "$opt\n    convergence<MolEnergyConvergence>: (";
    $opt = "$opt\n      cartesian = yes";
    $opt = "$opt\n      energy = \$..:..:mole";
    $opt = "$opt\n    )";
    $opt = "$opt\n  )\n";

    my $freq = "";
    if ($qcinput->frequencies()) {
        $freq = "% vibrational frequency input";
        $freq = "$freq\n  freq<MolecularFrequencies>: (";
        if ($symmetry ne "C1") {
            $freq = "$freq\n    point_group<PointGroup>: symmetry = $symmetry";
        }
        $freq = "$freq\n    molecule = \$:molecule";
        $freq = "$freq\n  )\n";
    }

    my $mpqcstart = sprintf ("mpqc: (\n  checkpoint = %s\n",
                             bool_to_yesno($qcinput->checkpoint()));
    $mpqcstart = sprintf ("%s  savestate = %s\n",
                          $mpqcstart,bool_to_yesno($qcinput->checkpoint()));
    $mpqcstart = sprintf ("%s  restart = %s\n",
                          $mpqcstart,bool_to_yesno($qcinput->restart()));
    if ($use_cints) {
        $mpqcstart = "$mpqcstart  integrals<IntegralCints>: ()\n";
    }
    my $mpqcstop = ")\n";
    my $emacs = "% Emacs should use -*- KeyVal -*- mode\n";
    my $warn = "% this file was automatically generated\n";
    my $lab = $qcinput->label();
    my $label = "";
    if (! $lab =~ /^\s*$/) {
        $label = "% label: $lab";
        $label =~ s/\n/\n% label: /g;
        $label = "$label\n";
    }
    "$emacs$warn$label$mol$basis$integrals$mpqcstart$coor$mole$opt$freq$mpqcstop";
}

sub mpqc_fixed_coor {
    my $self = shift;
    my $coorref = shift;
    my $result = "";
    $result = "\n    fixed<SetIntCoor>: [";
    while (scalar(@{$coorref}) != 0) {
        my $nextcoor = $self->mpqc_sum_coor("      ","",$coorref);
        $result = "$result\n$nextcoor";
    }
    $result = "$result\n    ]";
}

sub mpqc_followed_coor {
    my $self = shift;
    my $coorref = shift;
    sprintf "\n%s", $self->mpqc_sum_coor("    ","followed",$coorref);
}

sub mpqc_sum_coor {
    my $self = shift;
    my $space = shift;
    my $name = shift;
    my $coor = shift;
    my $result = "$space$name<SumIntCoor>:(";
    $result = "$result\n$space  coor: [";
    my @coef = ();
    do {
        $coef[$ncoor] = shift @{$coor};
        my $simple = $self->mpqc_coor($coor);
        $result = "$result\n$space    $simple";
        $ncoor = $ncoor + 1;
    } while($coor->[0] eq "+" && shift @{$coor} eq "+");
    $result = "$result\n$space  ]";
    $result = "$result\n$space  coef = [";
    my $i;
    foreach $i (0..$#coef) {
        $result = "$result $coef[$i]";
    }
    $result = "$result]";
    $result = "$result\n$space)";
    $result;
}

sub mpqc_coor {
    my $self = shift;
    my $coor = shift;
    my $type = shift @{$coor};
    if ($type eq "TORS") {
        return sprintf "<TorsSimpleCo>:(atoms = [%d %d %d %d])",
                       shift @{$coor},shift @{$coor},
                       shift @{$coor},shift @{$coor};
    }
    if ($type eq "BEND") {
        return sprintf "<BendSimpleCo>:(atoms = [%d %d %d])",
                       shift @{$coor},shift @{$coor},
                       shift @{$coor};
    }
    if ($type eq "STRE") {
        return sprintf "<StreSimpleCo>:(atoms = [%d %d])",
                        shift @{$coor},shift @{$coor};
    }
}

sub bool_to_yesno {
    if (shift) { return "yes"; }
    else { return "no"; }
}

sub append_reference {
    my $mole = shift;
    my $refmethod = shift;
    my $charge = shift;
    my $mult = shift;
    my $memory = shift;
    my $orthog_method = shift;
    my $lindep_tol = shift;
    my $docc = shift;
    my $socc = shift;
    my $guessbasis = shift;
    $mole = "$mole\n    reference<$refmethod>: (";
    $mole = "$mole\n      molecule = \$:molecule";
    $mole = "$mole\n      basis = \$:basis";
    $mole = "$mole\n      total_charge = $charge";
    $mole = "$mole\n      multiplicity = $mult";
    $mole = "$mole\n      memory = $memory";
    if ($orthog_method ne "" ) {
        $mole = "$mole\n      orthog_method = $orthog_method";
    }
    if ($lindep_tol ne "" ) {
        $mole = "$mole\n      lindep_tol = $lindep_tol";
    }
    if ($docc ne "") {$mole = "$mole\n      $docc";}
    if ($socc ne "") {$mole = "$mole\n      $socc";}
    if (! ($basis =~ /^STO/
           || $basis =~ /^MI/
           || $basis =~ /^\d-\d1G$/)) {
        $mole = "$mole\n      guess_wavefunction<$refmethod>: (";
        $mole = "$mole\n        molecule = \$:molecule";
        $mole = "$mole\n        total_charge = $charge";
        $mole = "$mole\n        multiplicity = $mult";
        if ($docc ne "") {$mole = "$mole\n        $docc";}
        if ($socc ne "") {$mole = "$mole\n        $socc";}
        $mole = "$mole\n        basis<GaussianBasisSet>: (";
        $mole = "$mole\n          molecule = \$:molecule";
        $mole = "$mole\n          name = \"$guessbasis\"";
        $mole = "$mole\n        )";
        $mole = "$mole\n        memory = $memory";
        $mole = "$mole\n      )";
    }
    $mole = "$mole\n    )";
    return $mole;
}

1;

