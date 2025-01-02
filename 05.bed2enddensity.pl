### Usage example
### intersectBed -a sample.bed -b region.interest.bed -wo | perl 05.bed2enddensity.pl

use strict;
use warnings;

if(@ARGV!=1)
{
	warn "Usage: <xxx.ol>\n";
	exit 1;
}


my %h;
my %TSS_end_5;
my %TSS_end_3;
my %TSS_cov;
my %TSS_len;
my $tot;
my $tot_cov;


my %P;
my %F;

my $win=80;

while(<>)
{
	# chr6    90525818        90525908        90      +       chr6    90524441        90534442        90
    # chr12   6776596 6776779 183     +       chr12   6767313 6777314 183

	chomp;
	my @F=split;
	my $cen=int(($F[9]+$F[10]+1)/2);
	my $dis2cen_5=$F[1]+1-$cen;
	my $dis2cen_3=$F[2]-$cen;
	
	my $L=$F[2]-$F[1];
	my $start= $dis2cen_5-$win;
	
#	if($F[13] eq "-")
#	{
#		($dis2cen_5, $dis2cen_3)=(-$dis2cen_3, -$dis2cen_5);
#		#$dis2cen_3=-$dis2cen_5;
#		$start= $dis2cen_3-$win;
#	}

	#for P
	for (1..($win-1))
	{
			$start++;
			$P{$start}++;
	}
	#for F
	for (1..($L-$win+1))
	{
			$start++;
			$F{$start}++;
	}
	#for P
	for (1..($win-1))
	{
			$start++;
			$P{$start}++;
	}	
	
	#unless(exists $h{"$F[6]:$F[7]:$F[8]"})
	{
		$TSS_end_5{$dis2cen_5}++;
		$TSS_end_3{$dis2cen_3}++;
		$TSS_cov{$dis2cen_5}++;
		$TSS_cov{$dis2cen_3}--;
		$TSS_len{int($dis2cen_5)}{short}++ if $F[3]<100;
		$TSS_len{int($dis2cen_5)}{long}++ if $F[3]>150;
		$TSS_len{int($dis2cen_5)}{mono}++ if $F[3]>166 && $F[3]<170;
		$TSS_len{int($dis2cen_3)}{short}++ if $F[3]<100;
		$TSS_len{int($dis2cen_3)}{long}++ if $F[3]>150;
		$TSS_len{int($dis2cen_3)}{mono}++ if $F[3]>166 && $F[3]<170;   
		$tot++;
		$tot_cov+=$F[-1]
	}
	#$h{"$F[6]:$F[7]:$F[8]"}++;
}

my $cov_tmp=0;

for my $pos (sort {$a<=>$b} keys %TSS_cov)
{
	$cov_tmp+=$TSS_cov{$pos};
	$TSS_cov{$pos}=$cov_tmp;
}

for my $pos (-250...250)
{
    my $a=$P{$pos} || 0; #partially covered
    my $b=$F{$pos} || 0; #fully covered
	
	print join("\t", $pos, $TSS_end_5{$pos}, $TSS_end_5{$pos}/$tot*1000, $TSS_end_3{$pos}, $TSS_end_3{$pos}/$tot*1000, $TSS_cov{$pos}/$tot_cov*1000, $TSS_len{$pos}{short}, ($TSS_len{$pos}{short}||0)/($TSS_len{$pos}{tot}||1), $TSS_len{$pos}{long}||0, ($TSS_len{$pos}{long}||0)/($TSS_len{$pos}{tot}||1), $a/($a+$b+1)*100), "\n";
}
