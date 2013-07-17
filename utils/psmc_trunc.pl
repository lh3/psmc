#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %opts = (n=>20, d=>0, s=>100);
getopts('n:d:s:', \%opts);

my $flag = 1;
my ($cut, $last_t, $first, @rs_lines);
while (<>) {
	$flag = ($1 == $opts{n}) if /^RD\s(\d+)/;
	next if $flag == 0;
	if (/^TR\s([\d\.]+)\s([\d\.]+)/) {
		$cut = $opts{d} * $opts{s} / $1;
		$last_t = 0;
		print;
	} elsif (/^RS/) {
		my @t = split;
		if ($last_t <= $cut && $cut < $t[2]) {
			$first = $#rs_lines;
		}
		push(@rs_lines, \@t);
		$last_t = $t[2];
	} elsif (/^PA\t(\S+)\s(.*)/) {
		my $rest = $2;
		my @time;
		# print RS lines
		for my $t (@rs_lines[$first .. $#rs_lines]) {
			$t->[1] -= $first;
			$t->[2] -= $cut;
			$t->[2] = $t->[2] > 0? $t->[2] : 0;
			$t->[2] = sprintf("%f", $t->[2]);
			push(@time, $t->[2]);
			print join("\t", @$t), "\n";
		}
		# generate the new pattern
		my @t = &parse_pattern($1);
		my $n_lambda = shift(@t);
		my $first_param = $t[$first];
		my ($last, $cnt, @s);
		$last = -1; $cnt = 0;
		for my $i ($first .. $#t) {
			if ($last != $t[$i]) {
				push(@s, $cnt) if $last >= 0;
				$last = $t[$i]; $cnt = 1;
			} else {
				++$cnt;
			}
		}
		push(@s, $cnt);
		# simplify the new pattern
		my $pat = '';
		$last = -1; $cnt = 0;
		for (@s) {
			if ($last != $_) {
				$pat .= $cnt > 1? "$cnt*$last+" : "$last+" if $last >= 0;
				$last = $_; $cnt = 1;
			} else {
				++$cnt;
			}
		}
		$pat .= $cnt > 1? "$cnt*$last" : "$last";
		# generate other parameters
		@t = split(/\s+/, $rest);
		my $theta = shift(@t);
		my $rho = shift(@t);
		shift(@t);
		print "PA\t$pat\t$theta\t$rho\t-1\t", join("\t", @t[$first_param..$#t]), "\t", join("\t", @time), "\n";
	} else { print; }
	$flag = 0 if /^\/\//;
}

sub parse_pattern {
	$_ = shift;
	s/\s+//g;
	@_ = split('\+');
	my $n_lambda = 0;
	my @stack;
	for (@_) {
		my ($x1, $x2) = (1, $_);
		$x1 = $1, $x2 = $2 if (/(\d+)\*(\d+)/);
		push(@stack, $x2) for (0 .. $x1-1);
		$n_lambda += $x1;
	}
	my @ret;
	for my $i (0 .. $#stack) {
		push(@ret, $i) for (0 .. $stack[$i]-1);
	}
	return ($n_lambda, @ret);
}
