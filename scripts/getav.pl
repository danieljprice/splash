#!/usr/bin/perl
# adds up the L2 errors from supersphplot output; calculates average
use strict;
use warnings;

use List::Util qw(sum);
my $file = ' ';
my $avg = 0;
if ($#ARGV < 0) {
   print "script which parses supersphplot output for L2 errors \n";
   print "and calculates the average. Written by D. Price. \n\n";
   die "Usage: $0 filename(s) \n";
}

foreach $file (@ARGV) {
   open my $fh, '<', $file or die "Can't open $file: $!";
   my @errors;

   while ( <$fh> ) {
      if ( my ($error) = m/L2 error\s+=\s+(\S*)\s+/ ) {
         #print "error = $error \n";
         push @errors, $error;
      }
   }

   my $nerrors = scalar(@errors);
   if ($nerrors > 0) {
      $avg = sum(@errors) / $nerrors;
   }
   else {
      $avg = 0;   
   }

   print "$file: Average of $nerrors errors: $avg\n";
}

#print "$_\n" for @errors;
