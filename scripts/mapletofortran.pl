#!/usr/bin/perl
#
# @(#) Perl script to extract and format Maple output to Fortran code
# @(#) that can be cut and pasted into a file
#
# @(#) (c) Daniel Price, Jan 2011, daniel.price@monash.edu
#
use Text::Wrap;

$Text::Wrap::columns=120;
$Text::Wrap::separator="&\n";

if ($#ARGV!=0) {
  die "Usage: $0 input.tex \n";
}
my ($infile) = @ARGV;
my @eqns = [];

open(INPUT, $infile) || die ("Can't open $infile $!");
my $ineq = 0;
my $equation = '';

while (<INPUT>) {
 #
 #--parse the maple .tex output, looking for
 #  equations enclosed as \mapleinline{}{}{}{ blah = \n ... }
 #  (only look for multi-line equations, ignore single line ones)
 #
 if ($ineq == 0) {
    if ( m/^\\mapleinline\{.*\}\{.*\}\{(.*=)\s*$/ ) {
       # start equation if beginning matches the above
       $equation=$1;
       $ineq = 1;
    }
 } elsif ( m/^(.*)\;\}.*$/ )  {
    # end the equation when we get to line containing } (close bracket)
    my $line = $1;
    $line =~ s/\s//g;
    $line =~ s/\^/\*\*/g;
    $line =~ s/\+/ \+ /g;
    $line =~ s/\-/ \- /g;

    # add the last line to the equation, editing as per below
    $equation = join("",$equation,$line);
    $ineq = 0;
    
    # wrap the final equation to 120 characters, with
    # Fortran continuation lines as above
    my $eqn = wrap('  ','        ',$equation);
    push @eqns, $eqn;
 } else {
    # append lines in the middle to the equation, stripping whitespace
    # converting ^ to ** and adding spaces around + and -
    my $line = $_;
    $line =~ s/\s//g;
    $line =~ s/\^/\*\*/g;
    $line =~ s/\+/ \+ /g;
    $line =~ s/\-/ \- /g;
    $equation = join("",$equation,$line);
  }
}

#
#--finally, print out the equations in reverse alphabetical order
#
sub backwards { $b cmp $a }

foreach $eqn (sort backwards @eqns) {
    print "$eqn \n\n";
}
print "\n\n";
