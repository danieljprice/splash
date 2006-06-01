#!/usr/bin/perl
# adds up the L2 errors from supersphplot output; calculates average
use strict;
use warnings;

use List::Util qw(sum);
my $file = ' ';

foreach $file (@ARGV) {
   open my $fh, '<', $file or die "Can't open $file: $!";
   my @errors;

   while ( <$fh> ) {
       my ($error) = m/L2 error\s+=\s+(\S*)\s+/;
       push @errors, $error;
   }

   my $avg = sum(@errors) / scalar(@errors);

   print "$file: Average error: $avg\n";
}

#print "$_\n" for @errors;
