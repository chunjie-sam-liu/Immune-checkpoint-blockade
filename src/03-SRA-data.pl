#! /usr/bin/perl -w
use strict;
use Getopt::Std;

use vars qw($opt_i $opt_o $opt_m);
getopts('i:o:m:');

open IN,"<$opt_i";
while(<IN>){
	chomp $_;
	$_ =~ /(...)(...).*/;
	print "lftp -e \"mirror /sra/sra-instant/reads/ByRun/sra/$1/$1$2/$_ /data/liucj/data/immune-checkpoint-blockade/ERP107734\" ftp-trace.ncbi.nih.gov>/data/liucj/data/immune-checkpoint-blockade/ERP107734/$_\n";
}
close IN;
