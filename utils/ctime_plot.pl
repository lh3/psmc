#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %opts = (d=>'ctime_plot', l=>10.0, h=>0.5, Y=>0, r=>'');
getopts('d:l:h:Y:r:L', \%opts);
die(qq(
Usage:   ctime_plot.pl [options] <in.ctime> [...]\n
Options: -d DIR      working directory [$opts{d}]
         -l FLOAT    length of the figure [$opts{l}]
         -h FLOAT    height of the figure [$opts{h}]
         -Y FLOAT    max Y-axis, 0 for auto [$opts{Y}]
         -r FILE     recombination map (.rho) [null]
\n)) if (-t STDIN && @ARGV == 0);
mkdir($opts{d}) unless (-d $opts{d});
my (%hash);
for my $f (0 .. $#ARGV) {
  my ($fh2, $last, $fh, @lc);
  $last = '';
  open($fh2, $ARGV[$f]);
  while (<$fh2>) {
	my @t = split;
	if ($t[0] ne $last) {
	  print {$fh} join("\t", @lc), "\n" if (defined $fh);
	  close($fh) if (defined $fh);
	  open($fh, ">$opts{d}/ct-$t[0].$f") || die;
	  $hash{$t[0]} = $t[2];
	  $last = $t[0];
	}
	$hash{$t[0]} = $t[2] if ($t[2] > $hash{$t[0]});
	printf {$fh} "%.6f\t$t[3]\n", $t[1]/1000000.0;
	@lc = ($t[2]/1000000.0, $t[3]);
  }
  print {$fh} join("\t", @lc), "\n";
  close($fh);
  close($fh2);
}

# load recombination map
my %hash_rho;
if ($opts{r}) {
  my ($fh2, $last, $fh);
  $last = '';
  $fh = undef;
  open($fh2, ($opts{r} =~ /\.gz$/)? "gzip -dc $opts{r} |" : $opts{r}) || die;
  while (<$fh2>) {
	my @t = split;
	if ($t[0] ne $last) {
	  close($fh) if (defined $fh);
	  open($fh, ">$opts{d}/rho-$t[0]") || die;
	  $hash_rho{$t[0]} = 1;
	  $last = $t[0];
	}
	printf $fh ("%.6f\t%.3f\n", $t[1]/1000000, $t[2]);
  }
  close($fh);
  close($fh2);
}

# plot
my $height = $opts{h};
my $yskip = 0.09;
for my $pre (sort keys %hash) {
  my $fh;
  open($fh, "| gnuplot") || die;
  my $h = $yskip + @ARGV * $height + ($hash_rho{$pre}? $height*0.7 : 0.0);
  my $maxy = $opts{Y}? qq(set yran [0:$opts{Y}];) : '';
  print $fh qq(set log y) if defined $opts{L};
  print $fh qq(
set t po eps co so 16;
set out "$opts{d}/$pre.eps";
set size $opts{l}, $h;
set xtics 1.0;
set mxtics 10;
set grid xtics;
set xran [0:*];
$maxy
set xlab "Coordinate (Mb)";
set lmar 12;
set rmar 2;
set multiplot;
);
  for my $f (0 .. $#ARGV) {
	my ($bmar, $tmar, $title);
	$title = $ARGV[$f];
	$title =~ s/\.ctime$//;
	$h = $f * $height + ($f == 0? 0.0 : $yskip);
	$bmar = ($f == 0)? 3 : 0;
	$tmar = ($f == $#ARGV && !$hash_rho{$pre})? 1 : 0;
	my $hei = $height + ($f == 0? $yskip : 0.0);
	my $str = '';
	$str .= qq(set xlab ""; set format x "";) if ($f > 0);
	print $fh qq(
set size $opts{l}, $hei;
set origin 0.0, $h;
set bmar $bmar;
set tmar $tmar;
set ylab "$title";
$str
plot "$opts{d}/ct-$pre.$f" t "" w st ls 1;
);
  }
  if ($hash_rho{$pre}) { # recombination map
	$h = @ARGV * $height + $yskip;
	my $tmph = $height * 0.7;
	print $fh qq(
set ylab "Recomb rate (cM/Mb)";
set size $opts{l}, $tmph;
set origin 0.0, $h;
set bmar 0;
set tmar 1;
set xlab ""; set format x "";
set yran [*:*];
plot "$opts{d}/rho-$pre" t "" w steps 3;
);
  }
  print $fh "exit;";
  close($fh);
}
