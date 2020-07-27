#! /usr/bin/perl
use PerlIO::gzip;
use Getopt::Long;
use File::Basename;
use FindBin '$Bin';
use Cwd qw(abs_path);
sub usage {
    print STDERR<< "USAGE";

Despriprion: BGI version's scIsoSeq gene expression count matrix generation program.
Author: shizhuoxing\@dingtalk.com
Data: 20200603
Usage: perl $0 -tgsbcumi bc_correction.xls -tmap filtered.tmap -topbc 2000 -sample test -outdir ./

Options:
	-ngsbc:			top rank of cellBC list in NGS data corresponding to the same library (optional)
	-tgsbcumi*:		tgs cellBC and umi correction result
	-tmap*:			tmap file from tgs FLNC to reference annotation gffcompare
	-topbc*:		top rank expression cellBC for tgs, mutex with -ngsbc
	-sample*:		sample name
	-outdir*:		output directory
	-help:			print this help
USAGE
    exit 0;
}

GetOptions(
	"ngsbc:s" => \$ngsbc,
	"tgsbcumi:s" => \$tgsbcumi,
	"tmap:s" => \$tmap,
	"topbc:s" => \$topbc,
	"sample:s" => \$sample,
	"outdir:s" => \$outdir,
	"help:s" => \$help
);
die &usage() if ((!defined $tgsbcumi) or (!defined $tmap) or (!defined $sample) or (defined $help));

if(defined $ngsbc){
	open IN,"$ngsbc";
	%ngs_bc=();
	while(<IN>){
		chomp;
		@a=();@a=split(/ /,$_);
		$ngs_bc{$a[0]}=$a[1];
	}
	close IN;
}

open IN,"$tgsbcumi";
%tgsbc=();%tgsumi=();
while(<IN>){
	chomp;
	@a=();@a=split(/\t/,$_);
	$tgsbc{$a[0]}=$a[1] if($a[3]>=0.95 and $a[5] ne "discarded" and $a[6] ne "discarded");
	$tgsumi{$a[0]}=$a[2] if($a[3]>=0.95 and $a[5] ne "discarded" and $a[6] ne "discarded");
}
close IN;

open IN,"$tmap";
$class_code="= c k m n j e o";
%matrix=();%tgsbccount=();
while(<IN>){
	chomp;
	@a=();@a=split(/\t/,$_);
	$a[3]=~s/.m1//;
	$bc="$tgsbc{$a[3]}-1";
	$umi="$tgsumi{$a[3]}";
	$matrix{$a[0]}{$bc}{$umi}++ if(defined $tgsbc{$a[3]} and $tgsumi{$a[3]} and $class_code=~/$a[2]/);
	$geneumi="$a[0].$umi";
	$tgsbccount{$bc}{$geneumi}++ if(defined $tgsbc{$a[3]} and $tgsumi{$a[3]} and $class_code=~/$a[2]/);
}
close IN;

if(!defined $ngsbc){
	%tgsbclist=();%tgs_bc=();
	@key1=();@key1=keys %tgsbccount;
	foreach $k1(@key1){
		@bc1=();@bc1=keys %{$tgsbccount{$k1}};
		$count1=@bc1;
		$tgsbclist{$k1}+=$count1;
	}
	
	$mark=0;
	foreach $k(sort {$tgsbclist{$b} <=> $tgsbclist{$a}} keys %tgsbclist){
		$mark++;
		$tgs_bc{$k}=1 if($mark<=$topbc);
		last if($mark>=$topbc);
	}
}

@key1=();@key1=keys %matrix;
@key2=();@key2=keys %ngs_bc if(defined $ngsbc);
@key2=keys %tgs_bc if(!defined $ngsbc);
open OUT1,">$outdir/$sample.gene.matrix";
foreach $k(@key2){
	$str.="\"$k\"\,";
}
$str=~s/\,$//;
print OUT1 "$str\n";

%tgsbc=();
foreach $k1(@key1){
	print OUT1 "\"$k1\"";
	foreach $k2(@key2){
		@key3=();@key3=keys %{$matrix{$k1}{$k2}};
		$count=@key3;
		$tgsbc{$k2}+=$count;
		print OUT1 ",$count";
	}
	print OUT1 "\n";
}
close OUT1;
