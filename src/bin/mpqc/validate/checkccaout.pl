#
eval 'exec perl $0 $*'
    if 0;

use Getopt::Long;

GetOptions("outputprefix=s" => \$outputprefix,
           "refprefix=s" => \$refprefix);

$all_passed=1;

foreach(@ARGV) {

  $file = $_;
  $file =~ s/(\.results)$//;

  # read in results
  $have_geom = $have_energy = $have_grad = 0;
  unless( open results, "$outputprefix/$file.results") {
    print "Couldn't open results for $outputprefix/$file.results\n";
  }

  $line = <results>;
  chomp($line);

  if( $line =~ /FINAL GEOMETRY\:/ ) {
    $have_geom=1;
    $line = <results>;
    if( $line =~ /^\n/ ) {
      print "problem with geometry\n";
      $have_geom=0; }
    while( !($line =~ /^\n/) ) {
      $line =~ /\s*(\S*)\s*(\S*)\s*(\S*)/;
      push(@geom,$1,$2,$3);
      $line = <results>;
    }
  }

  $line = <results>;
  if( $line =~ /FINAL ENERGY\:/ ) {
    $have_energy=1;
    $line = <results>;
    if( $line eq "\n" ) {
      print "problem with energy\n";
      $have_energy=0; }
    while( !($line =~ /^\n/) ) {
      $line =~ /\s*(\S*)/;
      $energy = $1;
      if( eof(results) ) {
        $line = "\n";
      }
      else {
        $line = <results>;
      }
    }
  }

  $line = <results>;
  if( $line =~ /FINAL GRADIENT\:/ ) {
    $have_grad=1;
    $line = <results>;
    if( $line =~ /^\n/ ) {
      print "problem with gradient\n";
      $have_grad=0; }
    while( !($line =~ /^\n/) ) {
      $line =~ /\s*(\S*)\s*(\S*)\s*(\S*)/;
      push(@grad,$1,$2,$3);
      if( eof(results) ) {
        $line = "\n";
      }
      else {
        $line = <results>;
      }
    }
  }
  close(results);

  # read in reference
  $have_ref_geom = $have_ref_energy = $have_ref_grad = 0;
  unless( open reference,
          "$refprefix/$file.results") {
    print "Couldn't open reference for $refprefix/$file.results\n";
  }
  $line = <reference>;
  chomp($line);

  if( $line =~ /FINAL GEOMETRY\:/ ) {
    $have_ref_geom = 1;
    $line = <reference>;
    if( $line =~ /^\n/ ) {
      print "problem with reference geometry\n";
      $have_ref_geom = 0; }
    while( !($line =~ /^\n/) ) {
      $line =~ /\s*(\S*)\s*(\S*)\s*(\S*)/;
      push(@ref_geom,$1,$2,$3);
      $line = <reference>;
    }
  }

  $line = <reference>;
  if( $line =~ /FINAL ENERGY\:/ ) {
    $have_ref_energy = 1;
    $line = <reference>;
    if( $line =~ /^\n/ ) {
      print "problem with reference energy\n";
      $have_ref_energy = 0; }
    while( !($line =~ /^\n/) ) {
      $line =~ /\s*(\S*)/;
      $ref_energy = $1;
      if( eof(reference) ) {
        $line = "\n";
      }
      else {
        $line = <reference>;
      }
    }
  }
  close(reference);


  # compare results to references
  open LOG, ">>$outputprefix/$file.diff" or
      print "Couldn't write to diff results for $file\n";
  $geom_passed = 1;
  if( $have_geom + $have_ref_geom != 2 ) {
    $geom_passed = 0; }
  @geom = reverse(@geom);
  @ref_geom = reverse(@ref_geom);
  print LOG "GEOMETRY DIFFERENCES WITH REFERENCE:\n";
  while(defined(@geom[0])) {
    $diff1 = pop(@geom) - pop(@ref_geom);
    $diff2 = pop(@geom) - pop(@ref_geom);
    $diff3 = pop(@geom) - pop(@ref_geom);
    printf(LOG "%20.9f%20.9f%20.9f\n", $diff1, $diff2, $diff3);
    if( abs($diff1) > $geom_tol |
        abs($diff2) > $geom_tol |
        abs($diff3) > $geom_tol   ) {
      $geom_passed=0;
    }
  }

  $energy_passed = 1;
  if( $have_energy + $have_ref_energy != 2 ) {
    $energy_passed = 0; }
  print LOG "\nENERGY DIFFERENCE WITH REFERENCE:\n";
  $diff1 = $energy - $ref_energy;
  printf(LOG "%20.9f\n", $diff1);
  if( abs($diff1) > $energy_tol ) {
    $energy_passed=0;
  }

  $grad_passed=1;
  if( $have_ref_grad ) {
    if( !$have_grad ) {
      $grad_passed = 0;
    }
    else {
      @grad = reverse(@grad);
      @ref_grad = reverse(@ref_grad);
      print LOG "\nGRADIENT DIFFERENCES WITH REFERENCE:\n";
      while(defined(@grad[0])) {
        $diff1 = pop(@grad) - pop(@ref_grad);
        $diff2 = pop(@grad) - pop(@ref_grad);
        $diff3 = pop(@grad) - pop(@ref_grad);
        printf(LOG "%20.9f%20.9f%20.9f\n", $diff1, $diff2, $diff3);
        if( abs($diff1) > $grad_tol |
            abs($diff2) > $grad_tol |
            abs($diff3) > $grad_tol   ) {
          $grad_passed=0;
        }
      }
    }
  }

  if( ($geom_passed + $energy_passed + $grad_passed) == 3 ) {
    push(@passed,"$file...passed");
  }
  else {
    push(@passed,"$file...failed");
    $all_passed = 0;
  }

  if(!$geom_passed) {
    print LOG
        "FAILURE: geometry not found or outside tolerance\n";
  }
  if(!$energy_passed) {
    print LOG
        "FAILURE: energy not found or outside tolerance\n";
  }
  if(!$grad_passed) {
    print LOG
        "FAILURE: gradient not found or outside tolerance\n";
  }
  close(LOG);
  @geom = @ref_geom = @grad = @ref_grad = ();
}

# display results
print "------------------------------------------------------------\n";
print "VERIFICATION RESULTS";
print "\n------------------------------------------------------------";

foreach $pass (@passed) {
  print "\n$pass";
}

if($all_passed) {
  print "\nALL TESTS PASSED\n";
}
else {
  print "\nFAILURE DURING VERIFICATION: examine results files for details\n";
}

# also output to file
open REP, ">./report.txt" or
    print "Couldn't write to report file\n";
print REP "------------------------------------------------------------\n";
print REP "VERIFICATION RESULTS";
print REP "\n------------------------------------------------------------";

foreach $pass (@passed) {
  print REP "\n$pass";
}

if($all_passed) {
  print REP "\nALL TESTS PASSED\n";
}
else {
  print REP
      "\nFAILURE DURING VERIFICATION: examine results files for details\n";
}
close(REP);

