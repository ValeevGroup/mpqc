#!/usr/bin/perl

# This script uses lwp-request from the libwww-perl package.
# The -o text option to lwp-request requires the HTML-Tree package
# and the HTML-Format package.

require 5.002;

$url = "http://www.emsl.pnl.gov/cgi-bin/ecce/basis_old.pl";
#$url = "http://localhost/cgi-bin/echo";

$SIG{'INT'} = 'dokill';

sub dokill {
  kill 9, $child if $child;
}

if ($#ARGV != 1) {
  printf "need two arguments: email_address and basis_name\n";
  exit;
}

$email = $ARGV[0];

$basis = "$ARGV[1]";

$basisreq = $basis;
$basisreq =~ s/\+/%2B/g;
$basisreq =~ s/\(/\%28/g;
$basisreq =~ s/\)/\%29/g;
$basisreq =~ s/,/\%2C/g;
$basisreq =~ s/ /\%20/g;

$atoms = "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr";
$atoms =~ s/ /+/g;

$email =~ s/@/%40/;

open(HOSTNAME,"hostname|");
$hostname = <HOSTNAME>;
close(HOSTNAME);
$hostname =~ s/\n//;

$data = sprintf "BasisSets=%s&Atoms=$atoms&Codes=NWChem&Optimize=on&ECP=on&Email=%s", $basisreq, $email;

$basisname = $basis;
$basisname =~ tr/A-Z/a-z/;
$basisname =~ tr/+/P/;
$basisname =~ tr/\*/S/;
$basisname =~ tr/\(/L/;
$basisname =~ tr/\)/R/;
$basisname =~ tr/,/_/;
$basisname =~ tr/ /_/;
$basisfile = "$basisname.nw";

printf "Basis will be put in %s.nw ...", $basisname;

open(HTTPD,"|lwp-request -m post -o text $url > $basisfile");
print HTTPD "$data";
close HTTPD;

printf " done.\n";
