#
eval 'exec perl $0 $*'
    if 0;

require QCParse;

##########################################################################

package Molecule;
$debug = 0;
$fltrx = "((?:-?\\d+|-?\\d+\\.\\d*|-?\\d*\\.\\d+)(?:[eE][+-]?\\d+)?)";

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

sub position {
    my $self = shift;
    my $i = shift;
    my $xyz = shift;
    $self->{"position"}[$i][$xyz];
}

1;
