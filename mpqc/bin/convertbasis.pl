#!/usr/bin/env perl

require AtomicBases;

my $basisname = shift;

my $bases = new AtomicBases;
$bases->set_name($basisname);

foreach my $file (@ARGV) {
    $bases->read_gaussian_file($file);
}

$bases->write_keyval(*STDOUT);


