#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %opts = (n=>20, s=>100, b=>0, e=>999);
getopts('n:s:b:e:', \%opts);

my ($flag, $DT) = (0, 0);
my ($theta, @t, @s);
while (<>) {
	$flag = 1 if (/^RD\t(\d+)/ && $1 == $opts{n});
	if ($flag) {
		$theta = $1 / $opts{s} if (/^TR\t(\S+)/);
		$DT = $1 * $theta if /^DT\t(\S+)/;
		($t[$1], $s[$1]) = ($2 * $theta, $3 * $theta) if /^RS\t(\S+)\t(\S+)\t(\S+)/;
		last if (/^\/\//);
	}
}
push(@t, 1000);
$opts{b} = $opts{b} > $DT? $opts{b} : 0;
$opts{e} = $opts{e} > $DT? $opts{e} : 0;
my $n = @s;
#for my $k (0 .. $#s) { print("$t[$k]\t$t[$k+1]\t$s[$k]\n"); }

my @alpha = (1);
for my $k (1 .. $n) {
	$alpha[$k] = $alpha[$k-1] * exp(-($t[$k] - $t[$k-1]) / $s[$k-1]);
}

my ($bk, $ek);
for my $k (1 .. $n) {
	$bk = $k - 1 if ($opts{b} < $t[$k] && $opts{b} >= $t[$k-1]);
	$ek = $k - 1 if ($opts{e} < $t[$k] && $opts{e} >= $t[$k-1]);
}

my $x = $alpha[$bk] * exp(-($opts{b} - $t[$bk]) / $s[$bk]) - $alpha[$ek] * exp(-($opts{e} - $t[$ek]) / $s[$ek]);

my $a1 = 0;
for (my $i = 0; $i < $bk; ++$i) {
	$a1 += $s[$i] * ($alpha[$i] - $alpha[$i+1]);
}
$a1 += $s[$bk] * $alpha[$bk] * (1 - exp(-($opts{b} - $t[$bk]) / $s[$bk]));

my $a2 = 0;
for (my $i = 0; $i < $ek; ++$i) {
	$a2 += $s[$i] * ($alpha[$i] - $alpha[$i+1]);
}
$a2 += $s[$ek] * $alpha[$ek] * (1 - exp(-($opts{e} - $t[$ek]) / $s[$ek]));

my $avg = ($a2 - $a1) / $x;

print "$avg\n";
