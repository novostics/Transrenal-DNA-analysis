use strict;
use warnings;

if(@ARGV!=1)
{
	warn "usage: <export.len.txt>\n";
	exit 1;
}
my %h;
my $tot;
while(<>)
{
	chomp;
	my @F=split;
	next if $F[3]>600;
	$h{$F[3]}++;
	$tot++;
}

print "#size\tcount\tpercentage\n";
for(0..600)
{
	$h{$_}=defined $h{$_} ? $h{$_} : 0;
	print join("\t", $_, $h{$_}, $h{$_}/$tot*100), "\n";
}
