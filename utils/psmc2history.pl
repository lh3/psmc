#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

&main;
exit;

sub main {
  my %opts = (n=>20, u=>-1);
  getopts('n:u:', \%opts);
  die("Usage: psmc2history.pl [-n $opts{n}] <in.psmc.par>\n") if (@ARGV == 0 && -t STDIN);
  my %h;
  $_ = <>;
  if (/^[A-Z][A-Z]/) {
  	my $flag = 0;
  	while (<>) {
		$flag = 1 if (/^RD.(\d+)/ and $1 == $opts{n});
		if ($flag and /^PA\t(.*)/) {
			$_ = $1;
			last;
		}
	}
	die unless $flag; 
  }
  &parse_param($_, \%h);
  my $a = $h{_};
  print "T $h{T}\nR $h{R}\n";
  print "N ", scalar(@{$h{_}}), "\n";
  if ($opts{u} < 0) {
	for my $x (0 .. @$a-1) {
		print "H $a->[$x][0]\t$a->[$x][1]\n";
	}
  } else {
  	my $N0 = $h{T} / 4 / $opts{u};
  	for my $x (0 .. @$a-1) {
		print "h ", $a->[$x][0] * 2 * $N0, "\t", $a->[$x][1] * $N0, "\n";
	}
  }
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
sub parse_param {
  my ($line, $ret) = @_;
  @_ = split(/\s+/, $line);
  $ret->{P} = shift; $ret->{T} = shift; $ret->{R} = shift;
  my ($n_lambda, @pat) = &parse_pattern($ret->{P});
  my $a = \@{$ret->{_}};
  $a->[$_][1] = shift for (0 .. $n_lambda-1);
  $a->[0][0] = 0;
  $a->[$pat[$_]+1][0] = $_[$_] for (0 .. @_-1);
  pop(@$a) if (!defined($a->[$n_lambda][1]));
}
