#! /usr/bin/perl
use PerlIO::gzip;
use Getopt::Long;
use File::Basename;
use FindBin '$Bin';
use Cwd qw(abs_path);
sub usage {
    print STDERR<< "USAGE";

Despriprion: BGI version's scIsoSeq isoform expression count matrix generation program.
Author: shizhuoxing\@dingtalk.com
Date: 20200603
Usage: perl $0 -tgsbcumi bc_correction.xls -group collapsed.group.txt -topbc 2000 -sample test -outdir ./

Options:
        -ngsbc:                 top rank of cellBC list in NGS data corresponding to the same library (optional)
        -tgsbcumi*:             tgs cellBC and umi correction result
        -group*:                collapsed.group.txt file from cDNA_cupcake 
        -topbc*:                top rank expression cellBC for tgs, mutex with -ngsbc
	-minUMIcount		filter out the minimum umi count isoform
        -sample*:               sample name
        -outdir*:               output directory
        -help:                  print this help
USAGE
    exit 0;
}

GetOptions(
        "ngsbc:s" => \$ngsbc,
        "tgsbcumi:s" => \$tgsbcumi,
        "group:s" => \$group,
        "topbc:s" => \$topbc,
        "minUMIcount:s" => \$minUMIcount,
        "sample:s" => \$sample,
        "outdir:s" => \$outdir,
        "help:s" => \$help
);
die &usage() if ((!defined $tgsbcumi) or (!defined $group) or (!defined $minUMIcount) or (!defined $sample) or (defined $help));

if(defined $ngsbc){
	open IN,"$ngsbc";
	%ngs_bc=();
	while(<IN>){
		chomp;
		@a=();@a=split(/ /,$_);
                $ngs_bc{$a[0]}=$a[1];
	}
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

open IN,"$group";
%matrix=();%tgsbccount=();
while(<IN>){
	chomp;
	@a=();@a=split;
	@b=();@b=split(/\,/,$a[1]);
	foreach $k(@b){
		$bc="$tgsbc{$k}-1";
		$umi=$tgsumi{$k};
		$matrix{$a[0]}{$bc}{$umi}++ if(defined $tgsbc{$k} and $tgsumi{$k});
		$isoumi="$a[0].$umi";
		$tgsbccount{$bc}{$isoumi}++ if(defined $tgsbc{$k} and $tgsumi{$k});
	}
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

open OUT1,">$outdir/$sample.isoform.matrix";
foreach $k2(@key2){
    $str.="\"$k2\",";
}
$str=~s/\,$//;
print OUT1 "$str\n";

foreach $k1(@key1){
    $exp_sum=0;
    $str="\"$k1\"";
    foreach $k2(@key2){
        if(defined $matrix{$k1}{$k2}){
		@key3=();@key3=keys %{$matrix{$k1}{$k2}};
            	$exp=@key3;
            	$str.="\,$exp";
            	$exp_sum+=$exp;
        }else{
            	$str.="\,0";
        }
    }
    print OUT1 "$str\n" if($exp_sum>=$minUMIcount);
}
