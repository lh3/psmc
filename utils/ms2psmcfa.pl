#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %opts = (s=>100);
getopts('s:', \%opts);
my ($k, $skip) = (0, $opts{s});
my @seq;

while (<>) {
  if (/^\@begin/) {
	++$k;
	$_ = <>;
	@seq = ();
	my $l = int($_/$skip) + 1;
	$seq[$_] = 0 for (0 .. $l-1);
  } elsif (/^(\d+)\s[01]+$/) {
	$seq[int($1/$skip)] = 1;
  } elsif (/\@end/) {
	print ">$k";
	for my $i (0 .. $#seq) {
	  print "\n" if ($i % 60 == 0);
	  print ($seq[$i]?'K':'T');
	}
	print "\n";
  }
}
