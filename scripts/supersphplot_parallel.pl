#!/bin/env perl
#
# farms out jobs to a list of machines, finding available CPU
#
use strict;
use warnings;
use IO::Handle;

my $run;
my $exe='~/ndspmhd/plot/ssupersphplot';
my $inputfile='input';
my $pgplotdev='none';
my $pgplotfile='pgplot';

if ($#ARGV < 1 ) {
   die "Usage: $0 nfilesperplot file1 [file2 file3 ... filen] \n also with a file called input \n";
}

my ($nfilesperplot,@files) = @ARGV;
my $pwd = `pwd`;
chomp($pwd);

if ($nfilesperplot < 1) { die "error: nfilesperplot must be > 0\n" };

my $nfiles=$#files + 1;
my $nruns=$nfiles/$nfilesperplot;
print "nfiles = $nfiles; nruns = $nruns\n";
#
#  get input from input file
#
my @inputs = `cat $inputfile` or die "error: input file does not exist\n";
#open my $fh,'<',$inputfile or die "error: input file does not exist\n";
my $line;
my $nline = 0;
my $devline = 0;
foreach $line (@inputs) {
    if ( $line =~ m|(/.+)\s| ) {
       $pgplotdev = $1;
       $devline = $nline;
    }
    $nline = $nline + 1;
}
#close $fh;

#
#--work out what device we are writing to and therefore what the filenames
#  should be
#
print "PGPLOT device = $pgplotdev $inputs[$devline]\n";
if ( $pgplotdev eq '/xw' or $pgplotdev eq '/aqt' ) { 
   die "must choose a non-interactive PGPLOT device\n";
}
my ($ext) = $pgplotdev =~ m|/(.+)|;
#print "extension = $ext \n";

#
#--farm out the jobs
#
my $filestart = 0;
my $fileend = $filestart + $nfilesperplot;
for ($run=1;$run<=$nruns;$run++) {
    my @inputstemp = @inputs;
    my $pgplotfile = "$files[$filestart].$ext";
    $inputstemp[$devline] = "$pgplotfile$pgplotdev\n";
    print "------------ run $run : $pgplotfile ------------\n";
    my @argsn=@files[$filestart..$fileend-1];
#    print "@inputstemp \n";
    open(INPUTF,"> input$run") || die("can't write temporary files");
    print INPUTF "@inputstemp \n";
    close(INPUTF);
#---here is the executable line---
    my $commandline = "cd $pwd; $exe @argsn < input$run \n";
    farmjob_xgrid( $commandline ) || die "error farming job \n";
#---------------------------------
    $filestart = $fileend;
    $fileend = $filestart + $nfilesperplot;
}
exit;

sub farmjob_xgrid {
    my $commandline=shift;
    print "command line is $commandline \n";
    my $xgridauth = "-hostname cytosine.ex.ac.uk -auth Kerberos";
    my $jobid = `xgrid $xgridauth -job submit $commandline` || die "xgrid not found \n";
    ($jobid) = $jobid =~ m/jobIdentifier\s+=\s+(\d+);/; # \s matches spaces (+ = at least one) \d decimals
    print "job id = $jobid \n";
    system "echo echo deleting... >> cleanup; echo xgrid $xgridauth -job delete -id $jobid >> cleanup";
}

sub farmjob_ssh {
    my (@machines) = `cat machinelist`;
    if ( $#machines < 1 ) {
       die "ERROR: no machines specified in file machinelist \n";
    } elsif ( $#machines < $nruns ) {
       print "WARNING: not enough machines given to run all jobs \n"
    }

    exit;
    my $machine;
    my $njobsrun = 1;
    my $jobsrun;
    my $ncpu = 1;
    my $njobspermachine = 3;
    my $n;

    # loop through all available machines looking for spare CPU
    foreach $machine (@machines) {
        chomp($machine);
        print "----------------- \n trying $machine \n";
    #    get load average for this machine using uptime command
        my $loadav = `ssh $machine uptime | cut -f5 -d':' | cut -f1 -d','`;
        chomp($loadav);
        print " load average last 1 minute = $loadav";
        if ( $loadav < 50.0 ) {
           print " ...OK \n";
           # run the job
           for ($n = 1;$n<=$njobspermachine;$n++){
    #           print "running $rootname$njobsrun on machine $machine at nice +19\n";
    #           system "ssh $machine 'cd $pwd/$rootname; nice +19 ./$ndim$SPMHD $rootname$njobsrun > $rootname$njobsrun.output &' ";
               $njobsrun = $njobsrun + 1;
               if ($njobsrun > $nruns) {
                  print "===========================================\n";
                  print "\n Hurrah! all jobs successfully submitted \n";
                  $jobsrun = $njobsrun - 1;
                  print " $jobsrun jobs run \n";
                  exit;
               }  
           } 
        } else {
           print " ...too busy \n";
        }
    }

    $jobsrun = $njobsrun - 1;
    print "=======================================================\n";
    print "WARNING: not enough machines available to run all jobs \n";
    print " $jobsrun jobs run \n";
    exit;
}

