#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

&main;
exit;

sub main {
  my %opts = (n=>2, L=>30000000, s=>100, u=>2.5e-8, R=>10, g=>25, d=>-1, r=>1); # u is useless actually
  getopts("n:L:s:u:r:R:Mg:d:r:", \%opts);
  die(qq/
Usage:   history2ms.pl [options] <in.psmc.par>

Options: -n INT    number of chromosome to simulate [$opts{n}]
         -L INT    length of each chromosome [$opts{L}]
         -s INT    skip used in psmc run [$opts{s}]
         -u FLOAT  neutral mutation rate [$opts{u}]
         -R FLOAT  recomb. rate in hotspots are FLOAT times larger [$opts{R}]
         -g INT    years per generation [$opts{g}]
         -d INT    divergence time [0]
         -r INT    # replicates [$opts{r}]
         -M        output macs command line
\n/) if (-t STDIN && @ARGV == 0);
  my ($theta, $rho, $n_lambda, $k, @rst, $N);
  while (<>) {
	if (/^T\s+(\S+)/) {
	  $theta = $1;
	  $N = $theta / $opts{s} / (4 * $opts{u});
	} elsif (/^R\s+(\S+)/) {
	  $rho = $1; $k = 0;
	} elsif (/^H\s+(\S+)\s+(\S+)/) {
	  $rst[$k]{B} = $1; $rst[$k++]{L} = $2;
	} elsif (/^h\s+(\S+)\s+(\S+)/) {
	  $rst[$k]{B} = $1 / 2 / $N / $opts{g}; $rst[$k++]{L} = $2 / $N;
	}
  }
  $n_lambda = $k;
# for my $x (0 .. $#rst) { print "$x\t$rst[$x]{B}\t$rst[$x]{L}\n"; }
  for my $p (@rst) {
	$p->{L} *= $N;
	$p->{B} *= 2 * $N;
  }
  # for Hudson's ms
  my $N0 = $rst[0]{L}; # present population size about 10000
  my $ms_theta = 4 * $N0 * $opts{u} * $opts{L};
  my $ms_rho = $ms_theta * ($rho / $theta);
  my ($macs_theta, $macs_rho) = ($ms_theta/$opts{L}, $ms_rho/$opts{L});
  my $ms_pop = '';
  my $n_tot = $opts{n};
  if ($opts{d} > 0) {
	  $opts{d} /= 4 * $N0;
	  $ms_pop = "-I 2 $opts{n} 1 -ej $opts{d} 2 1 -en 0 2 0.001";
	  ++$n_tot;
  }
  my $ms_cmd = defined($opts{M})? "macs $opts{n} $opts{L} -t $macs_theta -r $macs_rho" : "msHOT-lite $n_tot $opts{r} -t $ms_theta -r $ms_rho $opts{L} -l $ms_pop";
  for my $p (@rst) {
	$p->{B} /= 4 * $N0; $p->{L} /= $N0;
	next if ($p->{B} == 0.0);
	$ms_cmd .= sprintf(" -en %.4f 1 %.4f", $p->{B}, $p->{L});
  }
  print "$ms_cmd\n";
}
