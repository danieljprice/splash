#!/usr/bin/perl
# adds up the L2 errors from supersphplot output; calculates average
use strict;
use warnings;

use List::Util qw(sum);
my $file = ' ';
my $avg = 0;

foreach $file (@ARGV) {
   open my $fh, '<', $file or die "Can't open $file: $!";
   my @errors;

   while ( <$fh> ) {
      if (-m 'L2 error') {
         my ($error) = m/L2 error\s+=\s+(\S*)\s+/;
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
