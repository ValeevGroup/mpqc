#
eval 'exec perl $0 $*'
    if 0;

require AtomicBasis;

package AtomicBases;

$symbol_to_name{"H"}="hydrogen";
$symbol_to_name{"He"}="helium";
$symbol_to_name{"Li"}="lithium";
$symbol_to_name{"Be"}="beryllium";
$symbol_to_name{"B"}="boron";
$symbol_to_name{"C"}="carbon";
$symbol_to_name{"N"}="nitrogen";
$symbol_to_name{"O"}="oxygen";
$symbol_to_name{"F"}="fluorine";
$symbol_to_name{"Ne"}="neon";
$symbol_to_name{"Na"}="sodium";
$symbol_to_name{"Mg"}="magnesium";
$symbol_to_name{"Al"}="aluminum";
$symbol_to_name{"Si"}="silicon";
$symbol_to_name{"P"}="phosphorus";
$symbol_to_name{"S"}="sulfur";
$symbol_to_name{"Cl"}="chlorine";
$symbol_to_name{"Ar"}="argon";
$symbol_to_name{"K"}="potassium";
$symbol_to_name{"Ca"}="calcium";
$symbol_to_name{"Sc"}="scandium";
$symbol_to_name{"Ti"}="titanium";
$symbol_to_name{"V"}="vanadium";
$symbol_to_name{"Cr"}="chromium";
$symbol_to_name{"Mn"}="manganese";
$symbol_to_name{"Fe"}="iron";
$symbol_to_name{"Co"}="cobalt";
$symbol_to_name{"Ni"}="nickel";
$symbol_to_name{"Cu"}="copper";
$symbol_to_name{"Zn"}="zinc";
$symbol_to_name{"Ga"}="gallium";
$symbol_to_name{"Ge"}="germanium";
$symbol_to_name{"As"}="arsenic";
$symbol_to_name{"Se"}="selenium";
$symbol_to_name{"Br"}="bromine";
$symbol_to_name{"Kr"}="krypton";
$symbol_to_name{"Rb"}="rubidium";
$symbol_to_name{"Sr"}="strontium";
$symbol_to_name{"Y"}="yttrium";
$symbol_to_name{"Zr"}="zirconium";
$symbol_to_name{"Nb"}="niobium";
$symbol_to_name{"Mo"}="molybdenum";
$symbol_to_name{"Tc"}="technetium";
$symbol_to_name{"Ru"}="ruthenium";
$symbol_to_name{"Rh"}="rhodium";
$symbol_to_name{"Pd"}="palladium";
$symbol_to_name{"Ag"}="silver";
$symbol_to_name{"Cd"}="cadminium";
$symbol_to_name{"In"}="indium";
$symbol_to_name{"Sn"}="tin";
$symbol_to_name{"Sb"}="antimony";
$symbol_to_name{"Te"}="tellurium";
$symbol_to_name{"I"}="iodine";
$symbol_to_name{"Xe"}="xenon";
$symbol_to_name{"Cs"}="cesium";
$symbol_to_name{"Ba"}="barium";
$symbol_to_name{"La"}="lanthanium";
$symbol_to_name{"Ce"}="cerium";
$symbol_to_name{"Pr"}="praseodymium";
$symbol_to_name{"Nd"}="neodymium";
$symbol_to_name{"Pm"}="promethium";
$symbol_to_name{"Sm"}="samarium";
$symbol_to_name{"Eu"}="europium";
$symbol_to_name{"Gd"}="gadolinium";
$symbol_to_name{"Tb"}="terbium";
$symbol_to_name{"Dy"}="dysprosium";
$symbol_to_name{"Ho"}="holmium";
$symbol_to_name{"Er"}="erbium";
$symbol_to_name{"Tm"}="thulium";
$symbol_to_name{"Yb"}="ytterbium";
$symbol_to_name{"Lu"}="lutetium";
$symbol_to_name{"Hf"}="hafnium";
$symbol_to_name{"Ta"}="tantalum";
$symbol_to_name{"W"}="tungsten";
$symbol_to_name{"Re"}="rhenium";
$symbol_to_name{"Os"}="osmium";
$symbol_to_name{"Ir"}="iridium";
$symbol_to_name{"Pt"}="platinum";
$symbol_to_name{"Au"}="gold";
$symbol_to_name{"Hg"}="mercury";
$symbol_to_name{"Tl"}="thallium";
$symbol_to_name{"Pb"}="lead";
$symbol_to_name{"Bi"}="bismuth";
$symbol_to_name{"Po"}="polonium";
$symbol_to_name{"At"}="astatine";
$symbol_to_name{"Rn"}="radon";
$symbol_to_name{"Fr"}="francium";
$symbol_to_name{"Ra"}="radium";
$symbol_to_name{"Ac"}="actinium";
$symbol_to_name{"Th"}="thorium";
$symbol_to_name{"Pa"}="protactinium";
$symbol_to_name{"U"}="uranium";
$symbol_to_name{"Np"}="neptunium";
$symbol_to_name{"Pu"}="plutonium";
$symbol_to_name{"Am"}="americium";
$symbol_to_name{"Cm"}="curium";
$symbol_to_name{"Bk"}="berkelium";
$symbol_to_name{"Cf"}="californium";
$symbol_to_name{"Es"}="einsteinum";
$symbol_to_name{"Fm"}="fermium";
$symbol_to_name{"Md"}="mendelevium";
$symbol_to_name{"No"}="nobelium";
$symbol_to_name{"Lr"}="lawrencium";

@number_to_symbol = ();
$number_to_symbol[1] = "H";
$number_to_symbol[2] = "He";
$number_to_symbol[3] = "Li";
$number_to_symbol[4] = "Be";
$number_to_symbol[5] = "B";
$number_to_symbol[6] = "C";
$number_to_symbol[7] = "N";
$number_to_symbol[8] = "O";
$number_to_symbol[9] = "F";
$number_to_symbol[10] = "Ne";
$number_to_symbol[11] = "Na";
$number_to_symbol[12] = "Mg";
$number_to_symbol[13] = "Al";
$number_to_symbol[14] = "Si";
$number_to_symbol[15] = "P";
$number_to_symbol[16] = "S";
$number_to_symbol[17] = "Cl";
$number_to_symbol[18] = "Ar";
$number_to_symbol[19] = "K";
$number_to_symbol[20] = "Ca";
$number_to_symbol[21] = "Sc";
$number_to_symbol[22] = "Ti";
$number_to_symbol[23] = "V";
$number_to_symbol[24] = "Cr";
$number_to_symbol[25] = "Mn";
$number_to_symbol[26] = "Fe";
$number_to_symbol[27] = "Co";
$number_to_symbol[28] = "Ni";
$number_to_symbol[29] = "Cu";
$number_to_symbol[30] = "Zn";
$number_to_symbol[31] = "Ga";
$number_to_symbol[32] = "Ge";
$number_to_symbol[33] = "As";
$number_to_symbol[34] = "Se";
$number_to_symbol[35] = "Br";
$number_to_symbol[36] = "Kr";
$number_to_symbol[37] = "Rb";
$number_to_symbol[38] = "Sr";
$number_to_symbol[39] = "Y";
$number_to_symbol[40] = "Zr";
$number_to_symbol[41] = "Nb";
$number_to_symbol[42] = "Mo";
$number_to_symbol[43] = "Tc";
$number_to_symbol[44] = "Ru";
$number_to_symbol[45] = "Rh";
$number_to_symbol[46] = "Pd";
$number_to_symbol[47] = "Ag";
$number_to_symbol[48] = "Cd";
$number_to_symbol[49] = "In";
$number_to_symbol[50] = "Sn";
$number_to_symbol[51] = "Sb";
$number_to_symbol[52] = "Te";
$number_to_symbol[53] = "I";
$number_to_symbol[54] = "Xe";
$number_to_symbol[55] = "Cs";
$number_to_symbol[56] = "Ba";
$number_to_symbol[57] = "La";
$number_to_symbol[58] = "Ce";
$number_to_symbol[59] = "Pr";
$number_to_symbol[60] = "Nd";
$number_to_symbol[61] = "Pm";
$number_to_symbol[62] = "Sm";
$number_to_symbol[63] = "Eu";
$number_to_symbol[64] = "Gd";
$number_to_symbol[65] = "Tb";
$number_to_symbol[66] = "Dy";
$number_to_symbol[67] = "Ho";
$number_to_symbol[68] = "Er";
$number_to_symbol[69] = "Tm";
$number_to_symbol[70] = "Yb";
$number_to_symbol[71] = "Lu";
$number_to_symbol[72] = "Hf";
$number_to_symbol[73] = "Ta";
$number_to_symbol[74] = "W";
$number_to_symbol[75] = "Re";
$number_to_symbol[76] = "Os";
$number_to_symbol[77] = "Ir";
$number_to_symbol[78] = "Pt";
$number_to_symbol[79] = "Au";
$number_to_symbol[80] = "Hg";
$number_to_symbol[81] = "Tl";
$number_to_symbol[82] = "Pb";
$number_to_symbol[83] = "Bi";
$number_to_symbol[84] = "Po";
$number_to_symbol[85] = "At";
$number_to_symbol[86] = "Rn";
$number_to_symbol[87] = "Fr";
$number_to_symbol[88] = "Ra";
$number_to_symbol[89] = "Ac";
$number_to_symbol[90] = "Th";
$number_to_symbol[91] = "Pa";
$number_to_symbol[92] = "U";
$number_to_symbol[93] = "Np";
$number_to_symbol[94] = "Pu";
$number_to_symbol[95] = "Am";
$number_to_symbol[96] = "Cm";
$number_to_symbol[97] = "Bk";
$number_to_symbol[98] = "Cf";
$number_to_symbol[99] = "Es";
$number_to_symbol[100] = "Fm";
$number_to_symbol[101] = "Md";
$number_to_symbol[102] = "No";
$number_to_symbol[103] = "Lr";
$number_to_symbol[104] = "Rf";
$number_to_symbol[105] = "Ha";

%symbol_to_number = {};
foreach my $i (1..105) {
    $symbol_to_number{$number_to_symbol[$i]} = $i;
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

    $self->{"name"} = "";
    $self->{"atoms"} = {};
    $self->{"comment"} = {};
}

sub name {
    my $self = shift;
    return $self->{"name"};
}

sub set_name {
    my $self = shift;
    my $name = shift;
    $self->{"name"} = $name;
}

sub read_gaussian_file {
    my $self = shift;
    my $filename = shift;
    open(FILE,"<$filename");
    $self->read_gaussian(*FILE);
    close(FILE);
}


sub read_gaussian {
    my $self = shift;
    my $file = shift;
    while (<$file>) {
        if (/([A-Za-z]+)/) {
            my $sym = canonical_symbol($1);
            if (!exists($symbol_to_name{$sym})) {
                die "bad atomic symbol \"$sym\"";
            }
            $self->{"atoms"}->{$sym} = new AtomicBasis();
            $self->{"atoms"}->{$sym}->read_gaussian($file);
        }
    }
}

sub atoms {
    my $self = shift;
    return keys(%{$self->{"atoms"}});
}

sub comment {
    my $self = shift;
    return $self->{"comment"};
}

sub atom {
    my $self = shift;
    my $sym = canonical_symbol(shift);
    return $self->{"atoms"}->{$sym};
}

sub canonical_symbol {
    my $sym = shift;
    return uc(substr($sym,0,1)) . lc(substr($sym,1,length($sym)));
}

sub write_keyval_file {
    my $self = shift;
    my $filename = shift;
    open(FILE,">$filename");
    $self->write_keyval(*FILE);
    close(FILE);
};

sub write_keyval {
    my $self = shift;
    my $file = shift;
    my @atoms = $self->atoms();
    print $file "basis: (\n";
    foreach my $atom (sort { $symbol_to_number{$a} <=> $symbol_to_number{$b} }
                      @atoms) {
        printf $file " %s: \"%s\": [\n", $symbol_to_name{$atom}, $self->name();
        $self->atom($atom)->write_keyval($file);
        print $file " ]\n";
    }
    print $file ")\n";
}

1;
