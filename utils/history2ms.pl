#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

&main;
exit;

sub main {
  my %opts = (n=>2, L=>30000000, s=>100, u=>2.5e-8, r=>'', R=>10, g=>25); # u is useless actually
  getopts("n:L:s:u:r:R:Mg:", \%opts);
  die(qq/
Usage:   psmc2ms.pl [options] <in.psmc.par>

Options: -n INT    number of chromosome to simulate [$opts{n}]
         -L INT    length of each chromosome [$opts{L}]
         -s INT    skip used in psmc run [$opts{s}]
         -u FLOAT  neutral mutation rate [$opts{u}]
         -r FILE   recombination hotspot from HapMap [null]
         -R FLOAT  recomb. rate in hotspots are FLOAT times larger [$opts{R}]
         -g INT    years per generation [$opts{g}]
         -M        output macs command line
\n/) if (-t STDIN && @ARGV == 0);
  my $hs = &gen_hotspots(\%opts);
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
  my $ms_cmd = defined($opts{M})? "macs $opts{n} $opts{L} -t $macs_theta -r $macs_rho -T $hs" : "msHOT-lite $opts{n} 1 -t $ms_theta -r $ms_rho $opts{L} -T -l $hs";
  for my $p (@rst) {
	$p->{B} /= 4 * $N0; $p->{L} /= $N0;
	next if ($p->{B} == 0.0);
	$ms_cmd .= sprintf(" -en %.4f 1 %.4f", $p->{B}, $p->{L});
  }
  print "$ms_cmd\n";
}

sub gen_hotspots {
  my $opts = shift;
  return '' unless ($opts->{r});
  my %chr_len = (1=>247249719, 2=>242951149, 3=>199501827, 4=>191273063, 5=>180857866,
				 6=>170899992, 7=>158821424, 8=>146274826, 9=>140273252, 10=>135374737,
				 11=>134452384, 12=>132349534, 13=>114142980, 14=>106368585, 15=>100338915,
				 16=>88827254, 17=>78774742, 18=>76117153, 19=>63811651, 20=>62435964,
				 21=>46944323, 22=>49691432, X=>154913754);
  my ($chr, $start, $fh);
  my ($v, $n) = ('', 0);
  {
	my @t = keys(%chr_len);
	my $x = $t[int(rand()*@t)];
	$chr = "chr$x";
	$start = int(($chr_len{$x} - $opts->{L}) * rand());
  }
  $opts->{r} = "gzip -dc $opts->{r} |" if ($opts->{r} =~ /\.gz$/);
  open($fh, $opts->{r}) || die;
  my $last2 = 0;
  while (<$fh>) {
	my @t = split;
	next if ($t[0] ne $chr || $t[2] < $start || $t[3] >= $start + $opts->{L});
	next if ($last2 && $t[2] <= $last2);
	$last2 = $t[3];
	++$n;
	$v .= " " . join(" ", $t[2]-$start, $t[3]-$start, $opts->{R});
  }
  close($fh);
  return "-v $n$v";
}
