#!/usr/bin/perl -w

# Author: lh3

use strict;
use warnings;
use Getopt::Std;

my $version = "0.2.0";

my %opts = (u=>2.5e-8, 's'=>100, Y=>0, m=>5, X=>0, M=>'', x=>10000, n=>20, g=>25, f=>"Helvetica,16",
			w=>4, P=>"right top", T=>'', N=>0.0);
getopts('x:u:s:X:Y:RSGpm:n:M:N:g:f:w:P:T:L', \%opts);
die("
Usage:   psmc_plot.pl [options] <out.prefix> <in.psmc>\n
Options: -u FLOAT   absolute mutation rate per nucleotide [$opts{u}]
         -s INT     skip used in data preparation [$opts{s}]
         -X FLOAT   maximum generations, 0 for auto [0]
         -x FLOAT   minimum generations, 0 for auto [$opts{x}]
         -Y FLOAT   maximum popsize, 0 for auto [$opts{Y}]
         -m INT     minimum number of iteration [$opts{m}]
         -n INT     take n-th iteration (suppress GOF) [$opts{n}]
         -M titles  multiline mode [null]
         -f STR     font for title, labels and tics [$opts{f}]
         -g INT     number of years per generation [$opts{g}]
         -w INT     line width [$opts{w}]
         -P STR     position of the keys [$opts{P}]
         -T STR     figure title [null]
         -N FLOAT   false negative rate [0]
         -S         no scaling
         -L         show the last bin
         -p         convert to PDF (with epstopdf)
         -R         do not remove temporary files
         -G         plot grid
\n") if (@ARGV < 2);

my $prefix = shift(@ARGV);
my $scaling = defined($opts{S})? 0 : 1;
my (@data, $d, $N0, $skip, $Mseg, $Msize, $id, $min_ri, $do_store, $gof, $round, @FN, @nscale, @tscale, @alpha, $dt, @xshift, $N_scale, $fn0);

# initialize modifiers
if ($opts{M}) {
  my @t = split(/[,;]+/, $opts{M});
  for my $x (@t) {
	push(@FN, ($x =~ /=([^\s=*:@>]+)/)? $1 : 0.0);
	push(@nscale, ($x =~ /\*([^\s=*:@>]+)/? $1 : 1.0));
	push(@tscale, ($x =~ /:([^\s=*:@>]+)/? $1 : 0.0));
	push(@alpha, ($x =~ /\@([^\s=*:@>]+)/? $1 : 1.0));
	push(@xshift, ($x =~ />([^\s=*:@>]+)/? $1 : 0.0));
  }
  for (0 .. $#alpha) {
	$alpha[$_] = 2 * (2 + $alpha[$_]) / 3.0 / (1 + $alpha[$_]);
  }
}
$fn0 = $opts{N};

# load data

my $max_intv = 10000;
$id = $do_store = 0;
$gof = 'RI';
while (<>) {
  if (/^MM.*skip:(\d+)/) {
	$skip = $1 * $opts{s};
  } elsif (/^MM.*pattern:(\S+),/) {
  	my ($n_lambda, @pat) = &parse_pattern($1);
	for ($_ = $#pat - 1; $_ >= 0; --$_) {
		last if ($pat[$_] != $pat[$#pat]);
	}
	$max_intv = $_;
  } elsif (/^MM.*is_decoding:/) {
	$d = \%{$data[$id++]};
	$min_ri = 1e30; # reset
  } elsif (/^RD\s(\S+)/) {
	$Mseg = $Msize = $dt = 0;
	$round = $1;
  } elsif (/^(RI|GF)\s(\S+)/ && $1 eq $gof && $2 !~ /nan|inf/) {
	$do_store = 0;
	if ($round >= $opts{m} && $2 < $min_ri) {
	  $min_ri = $2;
	  $do_store = 1;
	}
	if ($opts{n} > 0) {
	  $do_store = ($round == $opts{n})? 1 : 0;
	}
  } elsif ($do_store && /^TR\s(\S+)\s(\S+)/) {
	($d->{T}, $d->{R}) = ($1/$skip, $2/$skip);
	$N_scale = 1.0;
	$N_scale /= 1.0 - $FN[$id-1] if (defined $FN[$id-1]);
	$N_scale /= 1.0 - $fn0;
	$N_scale /= $alpha[$id-1] if (defined $alpha[$id-1]);
	$N0 = $1/$skip / (4 * $opts{u}) * $N_scale;
	$d->{N0} = $N0;
	$d->{RI} = $min_ri;
  } elsif ($do_store && /^DT\s(\S+)/) {
  	$dt = $1;
  } elsif ($do_store && /^RS\s(\d+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)\s(\S+)/) { # psmc-0.6.0-5 or above
  	next if (!defined($opts{L}) && $1 > $max_intv);
	my $s = (defined $nscale[$id-1])? $nscale[$id-1] : 1.0;
	my $t = (defined $tscale[$id-1])? $tscale[$id-1] : 0.0;
	my $x = (defined $xshift[$id-1])? $xshift[$id-1] : 0.0;
	if ($scaling) {
		@{$d->{D}[$1]} = (2 * $N0 * ($2 + $dt) * (1.0-$t) * $opts{g} + $x, $3 * $N0 * $s / 10000, $4, $5, $6);
	} else {
		@{$d->{D}[$1]} = ($N_scale * $d->{T} * ($2 + $dt) * (1.0-$t) + $x, $N_scale * $3 * $d->{T} * $s * 1000, $4, $5, $6);
	}
	$Mseg = $4 if ($Mseg < $4);
	$Msize = 2 * $N0 * $2 if ($Msize < $2 * $N0);
  } elsif ($do_store && /^PA\s(.*)/) {
	$d->{PAR} = $1;
  } elsif ($do_store && /^\/\//) {
	$d->{Mseg} = $Mseg; $d->{Msize} = $Msize;
  }
}

# calculate tr_ratio and D_KL

my @misc;
$misc[0] = $data[0]{T} / $data[0]{R};
$misc[1] = $data[0]{RI};
my ($s1, $s2, $t1, $t2);
$s1 = $s2 = $t1 = $t2 = 0;
foreach my $i (0 .. @data-1) {
  $d = $data[$i];
  my $t = $d->{T} / $d->{R};
  $t1 += $t; $t2 += $t * $t;
  $s1 += $d->{RI}; $s2 += $d->{RI} * $d->{RI};
}
#$t1 /= @data; $t2 = sqrt(($t2 - @data * $t1 * $t1) / @data);
#$s1 /= @data; $s2 = sqrt(($s2 - @data * $s1 * $s1) / @data);
@misc[2..5] = ($t1, $s1, $t2, $s2);
# @misc = (tr, ri, avg_tr, avg_ri, dev_tr, dev_ri)

# write temporary file

my ($max_seg, $max_size, $fh);
$max_seg = $max_size = 0;
foreach my $i (0 .. @data-1) {
  $d = $data[$i];
  $max_seg = $d->{Mseg} if ($max_seg < $d->{Mseg});
  $max_size = $d->{Msize} if ($max_size < $d->{Msize});
  open($fh, ">$prefix.$i.txt") || die;
  foreach my $q (@{$d->{D}}) {
	print $fh join("\t", @$q), "\n";
  }
  close($fh);
}
# print .par
if ($opts{M}) {
  for (0 .. $#data) {
	my $d = $data[$_];
	open($fh, ">$prefix.$_.par") || die;
	print $fh "$d->{PAR}\n";
	close($fh);
  }
} else {
  open($fh, ">$prefix.par") || die;
  print $fh "$data[0]{PAR}\n";
  close($fh);
}

# plot
my $yran = ($opts{Y} > 0)? $opts{Y} : '*';
my $xran = ($opts{X} > 0)? $opts{X} : '*';
#my $title_str = sprintf('{/Symbol q} / {/Symbol r} = %.2f (%.2f +/- %.2f), D_{KL} = %.2e (%.2e +/- %.2e)',
#					$misc[0], $misc[2], $misc[4], $misc[1], $misc[3], $misc[5]);
my $keyconf = $opts{M}? "set key $opts{P}" : "set key off";
my $grid = $opts{G}? "set grid" : 'unset grid';
my $afont = qq/font "$opts{f}"/;
my $lw = qq/lw $opts{w}/;
my $ylab_aux = sprintf("%.1fx10^{-8}", $opts{u}/1e-8);

my ($xlab, $ylab);
if ($scaling) {
	$xlab = "Years (g=$opts{g}, {/Symbol m}=$ylab_aux)";
	$ylab = "Effective population size (x10^4)";
} else {
	$xlab = q[Time (scaled in units of 2{/Symbol m}T)];
	$ylab = q[Population size (scaled in units of 4{/Symbol m}N_e {/Symbol \264} 10^3)];
}

open($fh, "| tee $prefix.gp | gnuplot") || die;
print $fh qq(
  set size 1, 0.8;
  set xran [$opts{x}:$xran];
  set log x;
  set format x "10^{\%L}";
  set mxtics 10;
  set mytics 10;
  $grid;
  $keyconf;
  set xtics $afont;
  set ytics nomirror $afont;
  set xlab "$xlab" $afont;
  set t po eps enhance so co "Helvetica,16";
);
#print $fh qq(set title "$title_str";); # the title line
print $fh qq/set title "$opts{T}";/ if ($opts{T});
print $fh qq(
  set yran [0:$yran];
  set ylab "$ylab" $afont;
  set out "$prefix.eps";
  set style line 1 lt 1 lc rgb "#FF0000" $lw;
  set style line 2 lt 1 lc rgb "#00C000" $lw;
  set style line 3 lt 1 lc rgb "#0080FF" $lw;
  set style line 4 lt 1 lc rgb "#C000FF" $lw;
  set style line 5 lt 1 lc rgb "#00EEEE" $lw;
  set style line 6 lt 1 lc rgb "#C04000" $lw;
  set style line 7 lt 1 lc rgb "#C8C800" $lw;
  set style line 8 lt 1 lc rgb "#FF80FF" $lw;
  set style line 9 lt 1 lc rgb "#4E642E" $lw;
  set style line 10 lt 1 lc rgb "#800000" $lw;
  set style line 11 lt 1 lc rgb "#67B7F7" $lw;
  set style line 12 lt 1 lc rgb "#FFC127" $lw;
  plot );
if ($opts{M}) {
  my @titles = split(/[,;]/, $opts{M});
  for (0 .. $#titles) {
	$titles[$_] =~ s/=([^\s=*:@]+)/[$1]/;
	$titles[$_] =~ s/\@([^\s=*:@]+)/\(\{\/Symbol a\}=$1\)/;
  }
  foreach my $i (0 .. $#data) {
	print $fh qq("$prefix.$i.txt" u 1:2 t "$titles[$i]" w st ls ), $i + 1;
	print $fh ", " if ($i != $#data);
  }
} else {
  foreach my $i (1 .. @data-1) { print $fh qq("$prefix.$i.txt" u 1:2 w st not ls 2, ); }
  print $fh qq("$prefix.0.txt" u 1:2 t "popsize" w st ls 1);
  #print $fh qq(, "$prefix.0.txt" u 1:3 t "#segments" axis x1y2 w st 3);
  print $fh qq(;\n);
}
close($fh);

if (defined $opts{p}) {
  system("epstopdf $prefix.eps");
}

# remove files

unless (defined($opts{R})) {
  unlink <$prefix.*.txt>; unlink "$prefix.gp";
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
