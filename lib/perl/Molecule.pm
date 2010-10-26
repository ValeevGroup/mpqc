#
eval 'exec perl $0 $*'
    if 0;

use POSIX; # for acos
require QCParse;

##########################################################################

package Molecule;
$debug = 0;
$fltrx = "((?:-?\\d+|-?\\d+\\.\\d*|-?\\d*\\.\\d+)(?:[eE][+-]?\\d+)?)";

%z_to_sym = (
    "1" => "h",
    "2" => "he",
    "3" => "li",
    "4" => "be",
    "5" => "b",
    "6" => "c",
    "7" => "n",
    "8" => "o",
    "9" => "f",
    "10" => "ne",
    "11" => "na",
    "12" => "mg",
    "13" => "al",
    "14" => "si",
    "15" => "p",
    "16" => "s",
    "17" => "cl",
    "18" => "ar",
    "19" => "k",
    "20" => "ca",
    "21" => "sc",
    "22" => "ti",
    "23" => "v",
    "24" => "cr",
    "25" => "mn",
    "26" => "fe",
    "27" => "co",
    "28" => "ni",
    "29" => "cu",
    "30" => "zn",
    "31" => "ga",
    "32" => "ge",
    "33" => "as",
    "34" => "se",
    "35" => "br",
    "36" => "kr",
    "37" => "rb",
    "38" => "sr",
    "39" => "y",
    "40" => "zr",
    "41" => "nb",
    "42" => "mo",
    "43" => "tc",
    "44" => "ru",
    "45" => "rh",
    "46" => "pd",
    "47" => "ag",
    "48" => "cd",
    "49" => "in",
    "50" => "sn",
    "51" => "sb",
    "52" => "te",
    "53" => "i",
    "54" => "xe",
    "55" => "cs",
    "56" => "ba",
    "57" => "la",
    "58" => "ce",
    "59" => "pr",
    "60" => "nd",
    "61" => "pm",
    "62" => "sm",
    "63" => "eu",
    "64" => "gd",
    "65" => "tb",
    "66" => "dy",
    "67" => "ho",
    "68" => "er",
    "69" => "tm",
    "70" => "yb",
    "71" => "lu",
    "72" => "hf",
    "73" => "ta",
    "74" => "w",
    "75" => "re",
    "76" => "os",
    "77" => "ir",
    "78" => "pt",
    "79" => "au",
    "80" => "hg",
    "81" => "tl",
    "82" => "pb",
    "83" => "bi",
    "84" => "po",
    "85" => "at",
    "86" => "rn",
    "87" => "fr",
    "88" => "ra",
    "89" => "ac",
    "90" => "th",
    "91" => "pa",
    "92" => "u",
    "93" => "np",
    "94" => "pu",
    "95" => "am",
    "96" => "cm",
    "97" => "bk",
    "98" => "cf",
    "99" => "es",
    "100" => "fm",
    "101" => "md",
    "102" => "no",
    "103" => "lr",
    "104" => "rf",
    "105" => "ha"
);

%sym_to_z = (
    "h" => "1",
    "he" => "2",
    "li" => "3",
    "be" => "4",
    "b" => "5",
    "c" => "6",
    "n" => "7",
    "o" => "8",
    "f" => "9",
    "ne" => "10",
    "na" => "11",
    "mg" => "12",
    "al" => "13",
    "si" => "14",
    "p" => "15",
    "s" => "16",
    "cl" => "17",
    "ar" => "18",
    "k" => "19",
    "ca" => "20",
    "sc" => "21",
    "ti" => "22",
    "v" => "23",
    "cr" => "24",
    "mn" => "25",
    "fe" => "26",
    "co" => "27",
    "ni" => "28",
    "cu" => "29",
    "zn" => "30",
    "ga" => "31",
    "ge" => "32",
    "as" => "33",
    "se" => "34",
    "br" => "35",
    "kr" => "36",
    "rb" => "37",
    "sr" => "38",
    "y" => "39",
    "zr" => "40",
    "nb" => "41",
    "mo" => "42",
    "tc" => "43",
    "ru" => "44",
    "rh" => "45",
    "pd" => "46",
    "ag" => "47",
    "cd" => "48",
    "in" => "49",
    "sn" => "50",
    "sb" => "51",
    "te" => "52",
    "i" => "53",
    "xe" => "54",
    "cs" => "55",
    "ba" => "56",
    "la" => "57",
    "ce" => "58",
    "pr" => "59",
    "nd" => "60",
    "pm" => "61",
    "sm" => "62",
    "eu" => "63",
    "gd" => "64",
    "tb" => "65",
    "dy" => "66",
    "ho" => "67",
    "er" => "68",
    "tm" => "69",
    "yb" => "70",
    "lu" => "71",
    "hf" => "72",
    "ta" => "73",
    "w" => "74",
    "re" => "75",
    "os" => "76",
    "ir" => "77",
    "pt" => "78",
    "au" => "79",
    "hg" => "80",
    "tl" => "81",
    "pb" => "82",
    "bi" => "83",
    "po" => "84",
    "at" => "85",
    "rn" => "86",
    "fr" => "87",
    "ra" => "88",
    "ac" => "89",
    "th" => "90",
    "pa" => "91",
    "u" => "92",
    "np" => "93",
    "pu" => "94",
    "am" => "95",
    "cm" => "96",
    "bk" => "97",
    "cf" => "98",
    "es" => "99",
    "fm" => "100",
    "md" => "101",
    "no" => "102",
    "lr" => "103",
    "rf" => "104",
    "ha" => "105"
);

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
    my $mol = shift;

    $self->{"position"} = [];
    $self->{"element"} = [];
    my $i = 0;
    while ($mol =~ s/^\s*(\w+)\s+$fltrx\s+$fltrx\s+$fltrx\s*\n//) {
        my $sym = $1;
        my $x = $2;
        my $y = $3;
        my $z = $4;
        $self->{"element"}[$i] = $sym;
        $self->{"position"}[$i] = [ $x, $y, $z ];
        $i++;
    }

    $self->{"natom"} = $i;
}

sub string {
    my $self = shift;
    my $mol = "";
    for ($i = 0; $i < $self->n_atom(); $i++) {
        $mol = sprintf "%s  %s  %14.10f %14.10f %14.10f\n",
                       $mol, $self->element($i),
                       $self->position($i,0),
                       $self->position($i,1),
                       $self->position($i,2);
    }

    $mol;
}

sub n_atom {
    my $self = shift;
    printf "Molecule: returning natom = %d\n", $self->{"natom"} if ($debug);
    $self->{"natom"};
}

sub element {
    my $self = shift;
    my $i = shift;
    $self->{"element"}[$i];
}

sub z {
    my $self = shift;
    my $i = shift;
    $sym_to_z{lc($self->{"element"}[$i])};
}

sub position {
    my $self = shift;
    my $i = shift;
    my $xyz = shift;
    $self->{"position"}[$i][$xyz];
}

sub geometry {
    my $self = shift;
    my $geom = [];
    my $ngeom = 0;
    my $fence = $self->{"natom"}-1;
    foreach $i (0..$fence) {
        foreach $xyz (0..2) {
            $geom->[$ngeom] = $self->{"position"}[$i][$xyz];
            $ngeom = $ngeom + 1;
        }
    }
    $geom;
}

sub atom_xyz {
    my $self = shift;
    my $i = shift;
    my $geom = [];
    my $ngeom = 0;
    foreach $xyz (0..2) {
        $geom->[$ngeom] = $self->{"position"}[$i][$xyz];
        $ngeom = $ngeom + 1;
    }
    $geom;
}

sub dot {
    my $v1 = shift;
    my $v2 = shift;
    my $d = 0.0;
    foreach $xyz (0..2) {
        $d = $d + $v1->[$xyz] * $v2->[$xyz];
    }
    $d;
}

sub diff {
    my $v1 = shift;
    my $v2 = shift;
    my $diff = [];
    foreach $xyz (0..2) {
        $diff->[$xyz] = $v1->[$xyz] - $v2->[$xyz];
    }
    $diff;
}

sub vecstr {
    my $v = shift;
    sprintf "%12.8f %12.8f %12.8f", $v->[0], $v->[1], $v->[2];
}

# numbering starts at 1 for bond
sub bond {
    my $self = shift;
    my $a1 = shift() - 1;
    my $a2 = shift() - 1;
    my $d = diff($self->atom_xyz($a1),$self->atom_xyz($a2));
    #printf "v1 = %s\n", vecstr($self->atom_xyz($a1));
    #printf "v2 = %s\n", vecstr($self->atom_xyz($a2));
    #printf "diff = %s\n", vecstr($d);
    sqrt(dot($d,$d));
}

# numbering starts at 1 for bend
sub bend {
    my $self = shift;
    my $a1 = shift() - 1;
    my $a2 = shift() - 1;
    my $a3 = shift() - 1;
    my $diff12 = diff($self->atom_xyz($a1), $self->atom_xyz($a2));
    my $diff32 = diff($self->atom_xyz($a3), $self->atom_xyz($a2));
    POSIX::acos(dot($diff12,$diff32)/sqrt(dot($diff12,$diff12)*dot($diff32,$diff32))) * 180.0 / 3.14159265358979323846;
}

# numbering starts at 1 for tors
sub tors {
    my $self = shift;
    my $a1 = shift() - 1;
    my $a2 = shift() - 1;
    my $a3 = shift() - 1;
    my $a4 = shift() - 1;
    die "Molecule::tors not available";
}

1;
