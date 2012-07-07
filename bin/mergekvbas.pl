#!/usr/bin/perl

##############################################################################
#
#     BASIS SET FORMAT CONVERTER - merges 2 MPQC basis sets
#
#  Usage: mergekvbas.pl input_file_name1 input_file_name2 "new basis name" > output_file_name
#
##############################################################################
#

my $ifname1 = shift;
my $ifname2 = shift;
my $bsname  = shift;

open IFILE1, "<$ifname1" || die "could not open $ifname1";
open IFILE2, "<$ifname2" || die "could not open $ifname2";

my $token1;
while ( $token1 = <IFILE1> ) {
	last if $token1 =~ /basis\:/;
}

my $token2;
while ( $token2 = <IFILE2> ) {
	last if $token1 =~ /basis\:/;
}

printf STDOUT "$token1";

$token1 = <IFILE1>;
$token1 =~ s/ //g;
( $elemname1, $basisname1, $rest ) = split /\:/, $token1;
    
ELEM: while ($token2 = <IFILE2>) {
	
	$token2 =~ s/ //g;
	last ELEM if ($token2 =~ /\)/);
	( $elemname2, $basisname2, $rest ) = split /\:/, $token2;

	my $skip = ( $elemname1 ne $elemname2 );
    printf STDERR "skipped basis for $elemname2\n" if ($skip);
    printf STDERR "added basis for $elemname2\n" if (!$skip);

	if ( !$skip ) {
		printf STDOUT "  $elemname1: \"$bsname\": [\n";
		while ( $token1 = <IFILE1> ) {
			last if $token1 =~ / \]/;
			printf STDOUT "$token1";
		}
		
		$token1 = <IFILE1>;
        $token1 =~ s/ //g;
        ( $elemname1, $basisname1, $rest ) = split /\:/, $token1;
	}
	
	while ( $token2 = <IFILE2> ) {
		last if $token2 =~ / \]/;
		printf STDOUT "$token2" unless $skip;
	}
	printf STDOUT " ]\n" unless $skip;
}

printf STDOUT ")\n";

close IFILE1;
close IFILE2;
