#!/usr/bin/perl
#
# @(#) Perl script to extract release notes from src/splash.f90
#
my $start = 0;
my $out = "";
open(FILE,'src/splash.f90');
while (<FILE>) {
  my $line = $_;
  if ( m/(\d+\.\d+\.\d+)\s*:\s*\(\d+\/\d+\/\d+\)\s+(.*)$/ )  {
     if ($start == 1) {
        # print release notes as bullet points
        foreach (split /\; /, $out) {
           print "- $_\n";
        }
        close(FILE);
        exit();
     } else {
        $start = 1;
        $out = $2;
     }
  } elsif ($start == 1) {
     $line =~ m/\!\s+(.*)$/;  # match '!  some text' and get rid of ! and spaces
     $out = "${out} $1";      # append
  }
}
