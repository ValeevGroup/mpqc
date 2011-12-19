#!/usr/bin/perl

##############################################################################
#
#     BASIS SET FORMAT CONVERTER - GAUSSIAN_94 TO MPQC
#
###############################################################################
#
#  Used to convert a basis set formatted according to Gaussian 94 rules to
#  to the MPQC format. Reads g94 basis set from STDIN and prints out to STDOUT
#
#  Usage: g94_2_MPQC < input_file_name > output_file_name
#         g94_2_MPQC STDIN STDOUT
#
##############################################################################
# 

%elemname = (
  G  => ghost,
  H  => hydrogen,
  HE => helium,
  LI => lithium,
  BE => beryllium,
  B  => boron,
  C  => carbon,
  N  => nitrogen,
  O  => oxygen,
  F  => fluorine,
  NE => neon,
  NA => sodium,
  MG => magnesium,
  AL => aluminum,
  SI => silicon,
  P  => phosphorus,
  S  => sulfur,
  CL => chlorine,
  AR => argon,
  K  => POTASSIUM,
  CA => CALCIUM,
  SC => SCANDIUM,
  TI => TITANIUM,
  V  => VANADIUM,
  CR => CHROMIUM,
  MN => MANGANESE,
  FE => IRON,
  CO => COBALT,
  NI => NICKEL,
  CU => COPPER,
  ZN => ZINC,
  GA => GALLIUM,
  GE => GERMANIUM,
  AS => ARSENIC,
  SE => SELENIUM,
  BR => BROMINE,
  KR => KRYPTON,
  RB => RUBIDIUM,
  SR => STRONTIUM,
  Y  => YTTRIUM,
  ZR => ZIRCONIUM,
  NB => NIOBIUM,
  MO => MOLYBDENUM,
  TE => TECHNETIUM,
  RU => RUTHENIUM,
  RH => RHODIUM,
  PD => PALLADIUM,
  AG => SILVER,
  CD => CADMIUM,
  IN => INDIUM,
  SN => TIN,
  SB => ANTIMONY,
  TE => TELLURIUM,
  I  => IODINE,
  XE => XENON,
  CS => CESIUM,
  BA => BARIUM,
  LA => LANTHANUM,
  CE => CERIUM,
  PR => PRASEODYMIUM,
  ND => NEODYMIUM,
  PM => PROMETHIUM,
  SM => SAMARIUM,
  EU => EUROPIUM,
  GD => GADOLINIUM,
  TB => TERBIUM,
  DY => DYSPROSIUM,
  HO => HOLMIUM,
  ER => ERBIUM,
  TM => THULIUM,
  YB => YTTERBIUM,
  LU => LUTETIUM,
  HF => HAFNIUM,
  TA => TANATLUM,
  W  => TUNGSTEN,
  RE => RHENIUM,
  OS => OSMIUM,
  IR => IRIDIUM,
  PT => PLATINUM,
  AU => GOLD,
  HG => MERCURY,
  TL => THALLIUM,
  PB => LEAD,
  BI => BISMUTH,
  PO => POLONIUM,
  AT => ASTATINE,
  RN => RADON,
  FR => FRANCIUM,
  RA => RADIUM,
  U  => URANIUM,
  NP => NEPTUNIUM,
  PU => PLUTONIUM,
); 

# Find the first line containing BASIS keyword and parse basis name
BASIS: while (<STDIN>) {
    # Skip all lines until the first containing BASIS keyword
    ( !/BASIS|Basis|basis/) && next BASIS;
    ($scr, $basisname) = split "\"";
    last BASIS;
    }

my $pure = 0;
# make sure puream is 1 for correlation consistent and ano basis sets
# and 6-311g and sto-ng
if ($basisname =~ /cc-p/ || $basisname =~ /ANO/
    || $basisname =~ /^6-311G/
    || $basisname =~ /^6-311\+/
    || $basisname =~ /^STO-[1-9]G/
    ) {
  $pure = 1;
}


# Parse element name, contractions and primitives 
MAIN: while (<STDIN>) {
    # check if subsequent basis sets redefine basis name
    if (/BASIS|Basis|basis/) {
      ($scr, $basisname) = split "\"";
      next MAIN;
    }
    # Grab the element symbol and transform it to the standard
    # full element name
    chomp;
    @tokens = split " ";
    $elemsymb = $tokens[0];
    $elemsymb eq "" && last MAIN;
    $elemsymb =~ tr/[a-z]/[A-Z]/;
    print "basis:$elemname{$elemsymb}:\"$basisname\": [";

  # Loop over each contraction
  CONTR: while (<STDIN>) {
        /\*\*\*\*/ && last CONTR;
       ($am, $num_prim, $junk) = split " ";
       # Split SP contractions to S and P contractions 
       # Print in MPQC format
	if ($am =~ /SP/) {
          $num_prim2=$num_prim;
          # Read in each primitive
	  print "\n  (type: [am = S]\n    {exp coef:0} = {\n    ";
	  for($i=0; $i<$num_prim; $i++) {
             ($exp[$i], $s_coeff[$i], $p_coeff[$i] ) = split " ", <STDIN>;
             printf "  %30.12e   %30.12e", $exp[$i],$s_coeff[$i];
             $i != $num_prim - 1 && print "\n    ";
             } 
	  print "})"; # close the contraction
	  print "\n  (type: [(am = P puream = $pure)]\n    {exp coef:0} = {\n    ";
	  for($i=0; $i<$num_prim2; $i++) {
             printf "  %30.12e   %30.12e", $exp[$i],$p_coeff[$i];
             $i != $num_prim - 1 && print "\n    ";
	     } 
	  print "})"; # close the contraction
        }
        # Print the other contractions in MPQC format
        else{
          if ($am ne "S") {
            print "\n  (type: [(am = $am puream = $pure)]\n    {exp coef:0} = {\n    ";
          }
          else {
            print "\n  (type: [am = S]\n    {exp coef:0} = {\n    ";
          }
          # Read in each primitive
          for ($i=0; $i<$num_prim; $i++) {
             ($exp, $coeff) = split " ", <STDIN>;
             $exp =~ s/D/E/;
             $coeff =~ s/D/E/;
             printf "  %30.12e   %30.12e", $exp,$coeff;
             $i != $num_prim - 1 && print "\n    ";
          }    
          print "})"; # close the contraction
        }
  }
  print "\n]\n\n"; # close the basis set
}


