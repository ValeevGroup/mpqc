#!/usr/bin/perl

$port = 2080;
$them = 'www.emsl.pnl.gov';
$script = "/cgi-bin/run-bsform-post";
#$port = 80;
#$them = 'munin.ca.sandia.gov';
#$script = "/cgi-bin/echocat";

$AF_INET = 2;
$SOCK_STREAM = 1;

$SIG{'INT'} = 'dokill';

sub dokill {
  kill 9, $child if $child;
}

$sockaddr = 'S n a4 x8';


chop ($hostname = `hostname`);

($name, $aliases, $proto) = getprotobyname('tcp');
($name, $aliases, $port) = getservbyname($port,'tcp')
  unless $port =~ /^\d+$/;
($name, $aliases, $type, $len, $thataddr) = gethostbyname($them);
($name, $aliases, $type, $len, $thisaddr) = gethostbyname($hostname);

$this = pack($sockaddr, $AF_INET, 0, $thisaddr);
$that = pack($sockaddr, $AF_INET, $port, $thataddr);

if (socket(S, $AF_INET, $SOCK_STREAM, $proto)) {
  print "socket ok\n";
}
else {
  die $!;
}

if (bind(S, $this)) {
  print "bind ok\n";
}
else {
  die $!;
}

if (connect(S, $that)) {
  print "connect ok\n";
}
else {
  die $!;
}

select(S); $| = 1; select(STDOUT);

if ($#ARGV != 1) {
  printf "need two arguments: email_address and basis_name\n";
  exit;
}
$email = $ARGV[0];
$basis = $ARGV[1];
$basisurl = $basis;
$basisurl =~ s/\+/%2B/g;
$basisurl =~ s/\(/\%28/g;
$basisurl =~ s/\)/\%29/g;
$basisurl =~ s/,/\%2C/g;
printf "Basis is %s\n", $basis;
$data = sprintf "BasisSets=%s&Atoms=H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr&Optimize=on&ECP=on&Codes=NWCHEM&Email=%s", $basisurl, $email;
printf S "POST %s HTTP/1.0\n", $script;
printf S "Content-length: %d\n", length($data);
printf S "Content-type: application/x-www-form-urlencoded\n";
printf S "\n";
printf S "%s EOD\n", $data;

$name{"H"}="hydrogen";
$name{"He"}="helium";
$name{"Li"}="lithium";
$name{"Be"}="beryllium";
$name{"B"}="boron";
$name{"C"}="carbon";
$name{"N"}="nitrogen";
$name{"O"}="oxygen";
$name{"F"}="fluorine";
$name{"Ne"}="neon";
$name{"Na"}="sodium";
$name{"Mg"}="magnesium";
$name{"Al"}="aluminum";
$name{"Si"}="silicon";
$name{"P"}="phosphorus";
$name{"S"}="sulfur";
$name{"Cl"}="chlorine";
$name{"Ar"}="argon";
$name{"K"}="potassium";
$name{"Ca"}="calcium";
$name{"Sc"}="scandium";
$name{"Ti"}="titanium";
$name{"V"}="vanadium";
$name{"Cr"}="chromium";
$name{"Mn"}="manganese";
$name{"Fe"}="iron";
$name{"Co"}="cobalt";
$name{"Ni"}="nickel";
$name{"Cu"}="copper";
$name{"Zn"}="zinc";
$name{"Ga"}="gallium";
$name{"Ge"}="germanium";
$name{"As"}="arsenic";
$name{"Se"}="selenium";
$name{"Br"}="bromine";
$name{"Kr"}="krypton";
$name{"Rb"}="rubidium";
$name{"Sr"}="strontium";
$name{"Y"}="yttrium";
$name{"Zr"}="zirconium";
$name{"Nb"}="niobium";
$name{"Mo"}="molybdenum";
$name{"Tc"}="technetium";
$name{"Ru"}="ruthenium";
$name{"Rh"}="rhodium";
$name{"Pd"}="palladium";
$name{"Ag"}="silver";
$name{"Cd"}="cadminium";
$name{"In"}="indium";
$name{"Sn"}="tin";
$name{"Sb"}="antimony";
$name{"Te"}="tellurium";
$name{"I"}="iodine";
$name{"Xe"}="xenon";
$name{"Cs"}="cesium";
$name{"Ba"}="barium";
$name{"La"}="lanthanium";
$name{"Ce"}="cerium";
$name{"Pr"}="praseodymium";
$name{"Nd"}="neodymium";
$name{"Pm"}="promethium";
$name{"Sm"}="samarium";
$name{"Eu"}="europium";
$name{"Gd"}="gadolinium";
$name{"Tb"}="terbium";
$name{"Dy"}="dysprosium";
$name{"Ho"}="holmium";
$name{"Er"}="erbium";
$name{"Tm"}="thulium";
$name{"Yb"}="ytterbium";
$name{"Lu"}="lutetium";
$name{"Hf"}="hafnium";
$name{"Ta"}="tantalum";
$name{"W"}="tungsten";
$name{"Re"}="rhenium";
$name{"Os"}="osmium";
$name{"Ir"}="iridium";
$name{"Pt"}="platinum";
$name{"Au"}="gold";
$name{"Hg"}="mercury";
$name{"Tl"}="thallium";
$name{"Pb"}="lead";
$name{"Bi"}="bismuth";
$name{"Po"}="polonium";
$name{"At"}="astatine";
$name{"Rn"}="radon";
$name{"Fr"}="francium";
$name{"Ra"}="radium";
$name{"Ac"}="actinium";
$name{"Th"}="thorium";
$name{"Pa"}="protactinium";
$name{"U"}="uranium";
$name{"Np"}="neptunium";
$name{"Pu"}="plutonium";
$name{"Am"}="americium";
$name{"Cm"}="curium";
$name{"Bk"}="berkelium";
$name{"Cf"}="californium";
$name{"Es"}="einsteinum";
$name{"Fm"}="fermium";
$name{"Md"}="mendelevium";
$name{"No"}="nobelium";
$name{"Lr"}="lawrencium";


$basisname = $basis;
$basisname =~ tr/A-Z/a-z/;
$basisname =~ tr/+/P/;
$basisname =~ tr/\*/S/;
$basisname =~ tr/\(/L/;
$basisname =~ tr/\)/R/;
$basisname =~ tr/,/_/;
$basisname =~ tr/ /_/;
$basisfile = ">$basisname.kv";
#$basisfile = "|cat";
$atom = none;
$retrieve = 0;
$pure = 0;
open(BASIS, $basisfile);
printf BASIS "basis:(\n";
while (<S>) {
  #  print;
  #  next;
  GOTLINE:
  #printf "-----> %s\n", $_;
  if (/<pre>(BASIS.* +[^ ]* +)([A-Z]*)/) {
    $retrieve = 1;
    $line = "$1$2";
    if ($2 eq "SPHERICAL") {
      $pure = 1;
    }
    #printf "%s\n", $line;
    printf BASIS "%%%s\n", $line;
  }
  elsif (/^<p>END/) {
    $retrieve = 0;
  }
  elsif (/^<p>\#(.*)/) {
    printf BASIS "%%%s\n", $1;
  }
  elsif ($retrieve == 1) {
    /^<p>(.*)/;
    $_ = $1;
    #printf "%s\n", $_;
    if (/^ *([A-Z][a-z]*) +([A-Za-z]+)/) {
      if (!($1 eq $atom)) {
        if ($atom eq none) {
          $atom = $1;
          $am = $2;
          &start_atom;
          &start_shell;
        }
        else {
          &finish_shell;
          &finish_atom;
          $atom = $1;
          $am = $2;
          &start_atom;
          &start_shell;
        }
      }
      else {
        &finish_shell;
        $am = $2;
        &start_shell;
      }
    goto GOTLINE;
    }
    else {
      printf BASIS "   %s\n", $1;
    }
  }
}
if (!($atom eq none)) {
  &finish_shell;
  &finish_atom;
}
printf BASIS ")\n";
close(BASIS);
close(S);


sub start_atom {
  printf BASIS " %s: \"%s\": [\n", $name{$atom}, $basis;
}

sub finish_atom {
  printf BASIS " ]\n";
}

sub start_shell {
  while (<S>) {
    last;
  }
  @coefandexp = split(' ',$_);
  $ncoef = $#coefandexp - 1;
  ($amlower = $am) =~ tr/A-Z/a-z/;
  printf BASIS "  (type:";
  if ($amlower eq "sp") {
    printf BASIS " [am = p am = s]\n";
  }
  else {
    printf BASIS " [", $amlower;
    $icoef = 0;
    while ($icoef < $ncoef) {
      if ($icoef != 0) {
        printf BASIS " ";
      }
      if ($amlower eq "s" || $amlower eq "p" || !$pure) {
        printf BASIS "am = %s", $amlower;
      }
      else {
        printf BASIS "(am = %s puream = 1)", $amlower;
      }
      $icoef++;
    }
    printf BASIS "]\n", $amlower;
  }
  printf BASIS "   {exp";
  if ($amlower eq "sp") {
    printf BASIS " coef:1 coef:0";
  }
  else {
    $icoef = 0;
    while ($icoef < $ncoef) {
      printf BASIS " coef:%d", $icoef;
      $icoef++;
    }
  }
  printf BASIS "} = {\n";
}

sub finish_shell {
  printf BASIS "   })\n";
}
