#!/usr/bin/perl
#
# farms out jobs to a list of machines, finding available CPU
#
use strict;
use warnings;
use IO::Handle;
#---------------------------------------------------------
#  check args
#---------------------------------------------------------
if ($#ARGV < 1 ) {
   die "Usage: $0 nfilesperplot file1 [file2 file3 ... filen] \n also with a file called input \n";
}
#---------------------------------------------------------
#  set farming method options are ssh, xgrid
#---------------------------------------------------------
#my $farmusing='ssh';
my $farmusing='xgrid';
#---------------------------------------------------------
#  set name and location of the supersphplot executable
#---------------------------------------------------------
my $home=`cd; pwd -P`;
chomp ($home);
my $exe="$home/ndspmhd/plot/ssupersphplot";
my $pwd = `pwd`;
chomp($pwd);
#---------------------------------------------------------
#  set other default options
#---------------------------------------------------------
my $inputfile='input';
my $pgplotfile='pgplot';
my $pgplotdev='none'; # doesn't matter
my $pgplotdir= $ENV{'PGPLOT_DIR'}; # get this from the environment
my $ldpath= $ENV{'LD_LIBRARY_PATH'}; # get this from the environment
my $run;
#---------------------------------------------------------
#  get filenames and nfilesperplot from command line
#---------------------------------------------------------
my ($nfilesperplot,@files) = @ARGV;
#---------------------------------------------------------
#  work out how many files should be sent as command line arguments to each run
#---------------------------------------------------------
if ($nfilesperplot < 1) { die "error: nfilesperplot must be > 0\n" };
my $nfiles=$#files + 1;
my $nruns=$nfiles/$nfilesperplot;
print "nfiles = $nfiles; nruns = $nruns\n";
#---------------------------------------------------------
#  get generic input from input file
#---------------------------------------------------------
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

#---------------------------------------------------------
#  work out what device we are writing to and therefore 
#  what the filenames should be
#---------------------------------------------------------

print "PGPLOT device = $pgplotdev $inputs[$devline]\n";
if ( $pgplotdev eq '/xw' or $pgplotdev eq '/aqt' ) { 
   die "must choose a non-interactive PGPLOT device\n";
}
#---------------------------------------------------------
#  work out the filename extension from the device name
#---------------------------------------------------------
my ($ext) = $pgplotdev =~ m|/(.+)|;
#print "extension = $ext \n";

#---------------------------------------------------------
#  farm out the jobs
#---------------------------------------------------------
my $filestart = 0;
my $fileend = $filestart + $nfilesperplot;
my $n = 0;
for ($run=1;$run<=$nruns;$run++) {
    my @inputstemp = @inputs;
    my $num = sprintf("%04d",$run);
    my $pgplotfile = "pgplot_$num.$ext";
    #my $pgplotfile = "pgplot.$ext\_$run";
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
    print RUNSCR "source $home/.cshrc \n";
#    print RUNSCR "setenv LD_LIBRARY_PATH $ldpath \n";
    print RUNSCR "cd $pwd \n";
    print RUNSCR "$exe @argsn < input$run >& run$run.output \n";
    close(RUNSCR);
    system "chmod a+x run$run.csh";
    my $commandline = "./run$run.csh";
    if ( -s $pgplotfile ) {
       print "skipping $pgplotfile which already exists \n";
    } else {
      if ($farmusing eq 'ssh') {
         farmjob_ssh( $commandline, $n );    
      } elsif ($farmusing eq 'xgrid') {
         farmjob_xgrid( $commandline );
      } else {
         die "unknown farming method \n";
      }
    }
#---------------------------------
    $filestart = $fileend;
    $fileend = $filestart + $nfilesperplot;
}
exit;

#---------------------------------------------------------
#   subroutine to farm job via Apple's Xgrid
#---------------------------------------------------------

sub farmjob_xgrid {
    my $commandline=shift;
    my $xgridauth = "-hostname cytosine.ex.ac.uk -auth Kerberos -gid 0";
    my $jobid = `xgrid $xgridauth -job submit $commandline` || die "xgrid not found \n";
    ($jobid) = $jobid =~ m/jobIdentifier\s+=\s+(\d+);/; # \s matches spaces (+ = at least one) \d decimals
    print "farmed via xgrid: job id = $jobid \n";
    system "echo echo deleting... >> cleanxgrid; echo xgrid $xgridauth -job delete -id $jobid >> cleanxgrid";
    system "echo echo getting results... >> getresultsxgrid; echo xgrid $xgridauth -job results -id $jobid >> getresultsxgrid";
    sleep 1; # avoid multiple rapid-fire requests to xgrid server
}

#---------------------------------------------------------
#   subroutine to farm jobs via simple ssh commands
#   to a list of machines
#---------------------------------------------------------

sub farmjob_ssh {
    my ($commandline,$n)=@_;
#    my (@machines) = ('perky','tintin','homer','donald','butch','goofy','fred');
    my (@machines) = ('perky','tintin','obelix','homer','daphne','velma','taz','wallace','popeye','donald',
       'goofy','zippy','fred','butch');

    ##my (@machines) = `cat machinelist` or die "ERROR: for ssh version must list machines in file machinelist \n";
    my $nmachines = $#machines;
    if ( $nmachines < 1 ) {
       die "ERROR: no machines specified in file machinelist \n";
    }
    my $pwd = `pwd`;
    chomp($pwd);

    my $njobsrun = 0;
    my $loadavmax = 0.5;
    my $bs = '\\';

    # loop through all available machines looking for spare CPU
    while ($njobsrun < 1) {
        #--set which machine to try
        if ($n > $nmachines) { 
           print "waiting for machines to become free \n";
           sleep 30;
           $n = 0;
        };
        my $machine = $machines[$n];
        chomp($machine);
        
        my $loadav = 0.0;
        print "trying $machine: ";
        #--get load average for this machine using uptime command
        if ($loadavmax > 0.0) {
           my $uptime=`ssh $machine uptime`;
           ($loadav) = $uptime =~ m/load ave.+:\s(\d+\.\d+)./;
           print " load average last 1 minute = $loadav ";
        }
        #--set job running if not busy
        if ( $loadav < $loadavmax ) {
           print " ...OK \n";
           $njobsrun = $njobsrun + 1;
           # run the job
           print "running $commandline on machine $machine at nice +19\n";
           system "ssh -f $machine 'cd $pwd; nice +19 $commandline < /dev/null >& /dev/null &'";
        } else {
           print " ...too busy \n";
        }
        $n = $n + 1;
    }
    $_[1] = $n;
}
