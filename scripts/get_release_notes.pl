#!/usr/bin/perl
#
# @(#) Perl script to extract release notes from src/splash.f90
#
my $start = 0;
my $out = "";
my $do_all = 0;
if ( @ARGV[0] =~ m/all/ ) { $do_all = 1 };
open(FILE,'src/splash.f90');
while (<FILE>) {
  my $line = $_;
  if ( m/(\d+\.\d+\.\d+)\s*:.*$/)  {
     if ($start == 1) {
        # print release notes as bullet points
        foreach (split /\;/, $out) {
           print "-$_\n";
        }
        if ($do_all==1) {
           $out = "";
           print "\n:$1:\n\n";
        } else {
           close(FILE);
           exit()
        };
     } else {
        if ($do_all==1) {
           print "\n:$1:\n\n";
        }
        $start = 1;
        $out = $2;
     }
  } elsif ($start == 1) {
     # last entry, close on matching ---- line
     if (m/----/) {
        foreach (split /\;/, $out) {
           print "-$_\n";
        }
        close(FILE);
        exit()
     }
     $line =~ m/\!\s+(.*)$/;  # match '!  some text' and get rid of ! and spaces
     $out = "${out} $1";      # append
  }
}
