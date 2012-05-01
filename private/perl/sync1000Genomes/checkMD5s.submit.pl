#!/usr/bin/perl -w

use Getopt::Long;

sub usage {
    print "Usage: perl checkMD5s.submit.pl\n\t-chunk <chunk to check>\n";
    exit(1);
}


my $chunk = undef;
GetOptions( "chunk=s" => \$chunk);

usage() if ( !$chunk );

$cmd = "bsub -P checkMD5 -o $chunk.out -q gsa \"./checkMD5s.pl -ai $chunk\"";
system($cmd);
