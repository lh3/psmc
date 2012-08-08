#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %opts = (u=>2.5e-8, s=>100, N=>0.0);
getopts('u:s:rN:', \%opts); # 'r' for random generation
die("Usage: dec2ctime.pl [-u $opts{u}] [-s $opts{s}] [-N $opts{N}] <in.dec>\n") if (@ARGV == 0 && -t STDIN);

my $min_ri = 1e35;
my $do_store = 0;
my ($skip, $theta, $rho, $N0, @intv);
while (<>) {
  if (/^MM.*skip:(\d+)/) {
	$skip = $1 * $opts{s};
  } elsif (/^TR\s(\S+)\s(\S+)/) {
	$theta = $1 / $skip; $rho = $2 / $skip;
	$N0 = $theta / (4 * $opts{u}) / (1.0 - $opts{N});
  } elsif (/^RS\s+(\d+)\s+(\S+)/) {
	$intv[$1] = $2; $intv[@intv] = $intv[$#intv] + 3.0;
  } elsif (/^DC/) {
	my @t = split;
	$t[2] = ($t[2] - 1) * $skip + 1;
	$t[3] = $t[3] * $skip;
	$t[6] = defined($opts{r})? int(2 * $N0 * (rand($t[7] - $t[5]) + $t[5]))
	  : int(2 * $N0 * $t[6] + 0.5);
	print join("\t", @t[1..3,6]), "\n";
  }
}
