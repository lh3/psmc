#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %opts = (min_baseQ=>0);
getopts('q:', \%opts);
die("Usage: vcf2snp.pl [-q 0] <in.vcf>\n") if @ARGV < 1;

my %tab = (AC=>'M', AG=>'R', AT=>'W', CG=>'S', CT=>'Y', GT=>'K');

while (<>) {
	chomp;
	next if /^#/;
	my @t = split("\t");
	next unless $t[6] eq '.' || $t[6] eq 'PASS'; # skip if filtered
	next if $t[3] eq $t[4]; # skip if not variant
	my $v = uc($t[3] . $t[4]);
	next unless defined($tab{$v});
	print join("\t", $t[0], $t[1], $t[3], $tab{$v}), "\n";
}
