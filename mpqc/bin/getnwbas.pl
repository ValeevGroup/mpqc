#!/usr/bin/perl
require 5.002;
use Socket;

$port = 2080;
$remote = 'www.emsl.pnl.gov';
$script = "/cgi-bin/run-bsform-post";
#$port = 80;
#$remote = 'localhost';
#$script = "/cgi-bin/echocat";

$AF_INET = 2;
$SOCK_STREAM = 1;

$SIG{'INT'} = 'dokill';

sub dokill {
  kill 9, $child if $child;
}

if ($port =~ /\D/) { $port = getservbyname($port, 'tcp') }
die "No port" unless $port;
$iaddr   = inet_aton($remote)               || die "no host: $remote";

$paddr   = sockaddr_in($port, $iaddr);

$proto   = getprotobyname('tcp');
socket(S, PF_INET, SOCK_STREAM, $proto)  || die "socket: $!";
connect(S, $paddr)    || die "connect: $!";

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
$data = sprintf "BasisSets=%s&Atoms=H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr&Optimize=on&ECP=on&Codes=NWChem&Email=%s", $basisurl, $email;
printf S "POST %s HTTP/1.0\n", $script;
printf S "Content-length: %d\n", length($data);
printf S "Content-type: application/x-www-form-urlencoded\n";
printf S "\n";
printf S "%s EOD\n", $data;

$basisname = $basis;
$basisname =~ tr/A-Z/a-z/;
$basisname =~ tr/+/P/;
$basisname =~ tr/\*/S/;
$basisname =~ tr/\(/L/;
$basisname =~ tr/\)/R/;
$basisname =~ tr/,/_/;
$basisname =~ tr/ /_/;
$basisfile = ">$basisname.nw";
printf "Basis will be put in %s.nw ...", $basisname;
#$basisfile = "|cat";
open(BASIS, $basisfile);
while (<S>) {
    print BASIS;
}
close(BASIS);
close(S);
printf " done.\n";
