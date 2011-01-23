#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %opts = (s=>100, N=>0.0);
getopts('u:s:N:', \%opts);
die("Usage: decode2bed.pl [-s $opts{s}] [-N $opts{N}] <in.dec>\n") if (@ARGV == 0 && -t STDIN);

my ($skip, $theta, @intv);
while (<>) {
  if (/^MM.*skip:\s*(\d+)/) {
	$skip = $1 * $opts{s};
  } elsif (/^TR\s(\S+)\s(\S+)/) {
	$theta = ($1 / $skip) / (1. - $opts{N});
  } elsif (/^DC/) {
	my @t = split;
	$t[2] = ($t[2] - 1) * $skip;
	$t[3] = $t[3] * $skip;
	$t[5] = sprintf("%.8f", $t[5] * $theta);
	$t[6] = sprintf("%.8f", $t[6] * $theta);
	$t[7] = sprintf("%.8f", $t[7] * $theta);
	print join("\t", @t[1..5,7]), "\n";
  }
}
