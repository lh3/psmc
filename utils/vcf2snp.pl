#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %opts = (q=>0);
getopts('q:', \%opts);
die("Usage: vcf2snp.pl [-q 0] <in.vcf>\n") if @ARGV < 1;

my %tab = (AC=>'M', CA=>'M', AG=>'R', GA=>'R', AT=>'W', TA=>'W', CG=>'S', GC=>'S', CT=>'Y', TC=>'Y', GT=>'K', TG=>'K', AA=>'A', CC=>'C', GG=>'G', TT=>'T');

while (<>) {
	chomp;
	next if /^#/;
	my @t = split("\t");
	next unless $t[6] eq '.' || $t[6] eq 'PASS'; # skip if filtered
	next if $t[3] eq $t[4]; # skip if not variant
	next if length($t[3]) != 1; # not a SNP
	my @a = split(",", $t[4]);
	unshift(@a, $t[3]);
	my @s = split(":", $t[9]);
	@s = split(/[\/\|]/, $s[0]);
	my $v = uc($a[$s[0]] . $a[$s[1]]);
	next unless defined($tab{$v});
	print join("\t", $t[0], $t[1], $t[3], $tab{$v}), "\n";
}
