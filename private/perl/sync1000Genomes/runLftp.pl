#!/usr/bin/perl -w

# Runs lftp to pull down a file

use strict;
use Getopt::Long;

my $file = undef;
GetOptions( "file=s" => \$file);

if ( !$file ) {
    print "Usage: runLftp.pl\n\t-file \t<file>\n";
    exit(1);
}

chomp($file);

$file =~ m/ftp\/data\/(.*)\/(.*)\/.*/;

my $dir = "/humgen/1kg/DCC/ftp/data/$1/$2/";
mkdir($dir) unless(-d $dir);
chdir $dir;

open (TMP, '>lftp.txt');
print TMP "get ftp://ftp-trace.ncbi.nih.gov/1000genomes/$file\n";
print TMP "get ftp://ftp-trace.ncbi.nih.gov/1000genomes/$file.bai\n";
print TMP "get ftp://ftp-trace.ncbi.nih.gov/1000genomes/$file.bas\n";

# replace chrom11 with chrom20
$file =~ s/chrom11/chrom20/;

print TMP "get ftp://ftp-trace.ncbi.nih.gov/1000genomes/$file\n";
print TMP "get ftp://ftp-trace.ncbi.nih.gov/1000genomes/$file.bai\n";
print TMP "get ftp://ftp-trace.ncbi.nih.gov/1000genomes/$file.bas\n";
close (TMP);

my $cmd = "lftp -f lftp.txt";
system($cmd);

unlink "lftp.txt";
