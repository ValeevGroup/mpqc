#!/usr/bin/perl

$basisname = $ARGV[0];
$basisname =~ s/.nw$//;

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

$atom = none;
$retrieve = 0;
$pure = 0;
$pured = 0; # if $pure or $pured d's are pure
$puref = 1; # if $pure or $puref f's are pure; by default all f's are pure
# make sure puream is 1 for correlation consistent and ano basis sets
# and 6-311g and sto-ng
if ($basisname =~ /cc-p/ || $basisname =~ /ano/
    || $basisname =~ /^6-311g/
    || $basisname =~ /^sto-[1-9]g/
    ) {
    $pure = 1;
}
$basisname =~ tr/A-Z/a-z/;
$basisname =~ tr/+/P/;
$basisname =~ tr/\*/S/;
$basisname =~ tr/\(/L/;
$basisname =~ tr/\)/R/;
$basisname =~ tr/,/_/;
$basisname =~ tr/ /_/;
printf "Reading NWChem basis from %s.nw\n", $basisname;
printf "Writing MPQC basis to %s.kv\n", $basisname;
open(NWCHEMBASIS, "<$basisname.nw");
open(MPQCBASIS, ">$basisname.kv");
#open(MPQCBASIS, "|cat");
$firstatom=1;
$savedcomments="";
while (<NWCHEMBASIS>) {
  #  print;
  #  next;
  GOTLINE:
  #printf "-----> %s\n", $_;
  if (/<pre>(BASIS.* +[^ ]* +)([A-Z]*)/) {
    $retrieve = 1;
    $basis = $1;
    $line = "$1$2";
    $basis =~ s/^[^\"]*\"//; #"
    $basis =~ s/\"[^\"]*$//; #"
    printf "Basis = %s\n", $basis;
    if ($2 eq "SPHERICAL") {
      $pure = 1;
    }
    #printf "%s\n", $line;
    printf MPQCBASIS "%%%s\n", $line;
  }
  elsif (/^<p>END/) {
    $retrieve = 0;
  }
  elsif (/^<p>\#(.*)/) {
    my $comment = $1;
    $comment =~ s/ *$//;
    $savedcomments = sprintf("%s%%%s\n", $savedcomments, $comment);
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
      printf MPQCBASIS "   %s\n", $1;
    }
  }
}
if (!($atom eq none)) {
  &finish_shell;
  &finish_atom;
}
printf MPQCBASIS "%s", $savedcomments;
$savedcomments = "";
printf MPQCBASIS ")\n";
close(MPQCBASIS);
close(NWCHEMBASIS);


sub start_atom {
  if ($firstatom) {
    print MPQCBASIS "basis:(\n";
    $firstatom=0;
  }
  printf MPQCBASIS "%s", $savedcomments;
  $savedcomments = "";
  printf MPQCBASIS " %s: \"%s\": [\n", $name{$atom}, $basis;
}

sub finish_atom {
  printf MPQCBASIS " ]\n";
}

sub start_shell {
  printf MPQCBASIS "%s", $savedcomments;
  $savedcomments = "";
  while (<NWCHEMBASIS>) {
    last;
  }
  @coefandexp = split(' ',$_);
  $ncoef = $#coefandexp - 1;
  ($amlower = $am) =~ tr/A-Z/a-z/;
  printf MPQCBASIS "  (type:";
  if ($amlower eq "sp") {
    printf MPQCBASIS " [am = p am = s]\n";
  }
  else {
    printf MPQCBASIS " [", $amlower;
    $icoef = 0;
    while ($icoef < $ncoef) {
      if ($icoef != 0) {
        printf MPQCBASIS " ";
      }
      if ((($amlower eq "d") && $pured) || (($amlower eq "f") && $puref)) {
        printf MPQCBASIS "(am = %s puream = 1)", $amlower;
      }
      elsif ($amlower eq "s" || $amlower eq "p" || !$pure) {
        printf MPQCBASIS "am = %s", $amlower;
      }
      else {
        printf MPQCBASIS "(am = %s puream = 1)", $amlower;
      }
      $icoef++;
    }
    printf MPQCBASIS "]\n", $amlower;
  }
  printf MPQCBASIS "   {exp";
  if ($amlower eq "sp") {
    printf MPQCBASIS " coef:1 coef:0";
  }
  else {
    $icoef = 0;
    while ($icoef < $ncoef) {
      printf MPQCBASIS " coef:%d", $icoef;
      $icoef++;
    }
  }
  printf MPQCBASIS "} = {\n";
}

sub finish_shell {
  printf MPQCBASIS "   })\n";
}
