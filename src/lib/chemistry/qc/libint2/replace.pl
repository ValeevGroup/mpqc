#!/usr/bin/perl

my $what = "Cints";
my $with = "Libint2";

my $files = "Makefile\n";
$files .= `ls *.cc`;
$files .= `ls *.h`;

$_ = $files;
my @files = split;

printf STDOUT "@files\n";

foreach my $file (@files) {
  system("sed s/$what/$with/g $file > crap.txt");
  system("mv crap.txt $file");
}

