#!/usr/bin/perl
#
# farms out jobs to a list of machines, finding available CPU
#
use strict;
use warnings;
use IO::Handle;

my $run;
my $home=`cd; pwd -P`;
chomp ($home);
my $exe="$home/ndspmhd/plot/ssupersphplot";
my $inputfile='input';
my $pgplotdev='none';
my $pgplotfile='pgplot';
my $pgplotdir= $ENV{'PGPLOT_DIR'}; # get this from the environment

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
    my $num = sprintf("%05d",$run);
    my $pgplotfile = "pgplot_$num.$ext";
    $inputstemp[$devline] = "$pgplotfile$pgplotdev\n";
    print "------------ run $run : $pgplotfile ------------\n";
    my @argsn=@files[$filestart..$fileend-1];
#---write the input file
    open(INPUTF,"> input$run") || die("can't write temporary files");
    foreach $line (@inputstemp) {
       chomp($line);
       print INPUTF "$line \n";
    }
    close(INPUTF);
#---write the executable script---   
    open(RUNSCR,"> run$run.csh") || die("can't write run script");
    print RUNSCR "#!/bin/tcsh \n";
    print RUNSCR "setenv PGPLOT_DIR $pgplotdir \n";
    print RUNSCR "cd $pwd \n";
    print RUNSCR "$exe @argsn < input$run >& run$run.output \n";
    close(RUNSCR);
    system "chmod a+x run$run.csh";
    my $commandline = "./run$run.csh";
    farmjob_xgrid( $commandline );
#---------------------------------
    $filestart = $fileend;
    $fileend = $filestart + $nfilesperplot;
}
exit;

sub farmjob_xgrid {
    my $commandline=shift;
    my $xgridauth = "-hostname cytosine.ex.ac.uk -auth Kerberos";
    my $jobid = `xgrid $xgridauth -job submit $commandline` || die "xgrid not found \n";
    ($jobid) = $jobid =~ m/jobIdentifier\s+=\s+(\d+);/; # \s matches spaces (+ = at least one) \d decimals
    print "farmed via xgrid: job id = $jobid \n";
    system "echo echo deleting... >> cleanxgrid; echo xgrid $xgridauth -job delete -id $jobid >> cleanxgrid";
    system "echo echo getting results... >> getresultsxgrid; echo xgrid $xgridauth -job results -id $jobid >> getresultsxgrid";
    sleep 1; # avoid multiple rapid-fire requests to xgrid server
}

sub farmjob_ssh {
    my ($commandline,$runs)=@_;
    my (@machines) = `cat machinelist` or die "ERROR: for ssh version must list machines in file machinelist \n";
    my $nmachines = $#machines;
    if ( $nmachines < 1 ) {
       die "ERROR: no machines specified in file machinelist \n";
    }

    my $machine;
    my $njobsrun = 0;
    my $loadavmax = 0.5;
    my $n;

    # loop through all available machines looking for spare CPU
    while ($njobsrun < $runs) {
        
        chomp($machine);
        my $loadav = 0.;
        print "----------------- \n trying $machine \n";
    #    get load average for this machine using uptime command
        if ($loadavmax > 0.0) {
           my $uptime=`ssh $machine uptime`;
           ($loadav) = $uptime =~ m/load averages:(\d+\.\d+)\s/;
           print " load average last 1 minute = $loadav";
        }
        if ( $loadav < $loadavmax ) {
           print " ...OK \n";
           $njobsrun = $njobsrun + 1;
           # run the job
    #           print "running $rootname$njobsrun on machine $machine at nice +19\n";
    #           system "ssh $machine 'cd $pwd/$rootname; nice +19 ./$ndim$SPMHD $rootname$njobsrun > $rootname$njobsrun.output &' ";
        } else {
           print " ...too busy \n";
        }

    }
    print "===========================================\n";
    print "\n Hurrah! all jobs successfully submitted \n";
    print " $njobsrun jobs run \n";
}
