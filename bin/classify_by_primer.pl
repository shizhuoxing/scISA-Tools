#! /usr/bin/perl
use Getopt::Long;
use File::Basename;
use FindBin '$Bin';
use Cwd qw(abs_path);
sub usage {
    print STDERR<< "USAGE";

Despriprion: BGI version's full-length transcript detection algorithm for BGI patented linking scIsoSeq library construction protocol.
Author: shizhuoxing\@dingtalk.com
Date: 20191231
Usage: classify_by_primer -blastm7 mapped.m7 -ccsfq ccs.fq -min_primerlen 16 -min_isolen 50 -outdir ./

Options:
	-blastm7*:		result of primer blast to ccs.fa in blast -outfmt 7 format
	-ccsfq*:		the ccs.fq you want to classify to get full-length transcript
	-min_primerlen*:	the minimum primer alignment length in ccs.fa
	-min_seqlen*:		the minimum output's transcript length
	-outdir*:		output directory
	-help:			print this help
USAGE
    exit 0;
}

GetOptions(
	"blastm7:s" => \$blastm7,
	"ccsfq:s" => \$ccsfq,
	"min_primerlen:s" => \$min_primerlen,
	"min_seqlen:s" => \$min_seqlen,
	"outdir:s" => \$outdir,
	"help:s" => \$help
);
die &usage() if ((!defined $blastm7) or (!defined $ccsfq) or (!defined $min_primerlen) or (!defined $min_seqlen) or (!defined $outdir) or (defined $help));

open IN,"$blastm7";
%m7=();
while(<IN>){
	chomp;
	@a=();@a=split(/\t/,$_);$slen=@a;
	$m7{$a[0]}{$a[6]}="$a[1]\t$a[2]\t$a[3]\t$a[7]\t-" if($a[1] eq "primer_F" and $a[8]>$a[9] and $a[8]>=24 and $slen==12 and $a[3]>=$min_primerlen);
	$m7{$a[0]}{$a[6]}="$a[1]\t$a[2]\t$a[3]\t$a[7]\t+" if($a[1] eq "primer_F" and $a[8]<$a[9] and $a[9]>=24 and $slen==12 and $a[3]>=$min_primerlen);
	$m7{$a[0]}{$a[6]}="$a[1]\t$a[2]\t$a[3]\t$a[7]\t-" if($a[1] eq "primer_S" and $a[8]>$a[9] and $a[8]>=18 and $slen==12 and $a[3]>=$min_primerlen);
	$m7{$a[0]}{$a[6]}="$a[1]\t$a[2]\t$a[3]\t$a[7]\t+" if($a[1] eq "primer_S" and $a[8]<$a[9] and $a[9]>=18 and $slen==12 and $a[3]>=$min_primerlen);
}
close IN;

open IN,"$ccsfq";
open FL1,">$outdir/isoseq_flnc.BarcodeUMI.fastq";
open FL2,">$outdir/isoseq_flnc.Transcript.fastq";
open BARCODE,">$outdir/isoseq.PrimerStat.csv";
print BARCODE "SeqID\,SeqClassify\,SeqLength\,Primer1\,Primer1_BlastIdentity\,Primer1_BlastLength\,Primer2\,Primer2_BlastIdentity\,Primer2_BlastLength\,cellBC\,UMI\,PolyAtailLength\,PolyAtail\n";
while($id=<IN>){
	chomp $id;
	@a=();@a=split(/ /,$id);
	$name=$a[0];
	$name=~s/^@//;
	$tname=$name;$tname=~s/\/ccs//;

	$ccsseq=<IN>;chomp $ccsseq;
	$mk=<IN>;chomp $mk;
	$ccsqv=<IN>;chomp $ccsqv;

	$ccslen=length($ccsseq);
	@key1=();@key1=keys %{$m7{$name}};
	@key2=();@key2=sort{$a<=>$b}@key1;
	%sign=();%p=();$mm=0;$mx=@key2;
	foreach $k(@key2){
		$sign{$mm}=0;@{$p{$mm}}=();@{$p{$mm}}=split(/\t/,$m7{$name}{$k});
		
		if($mm>0){

			if($p{$mm-1}[0] eq "primer_S" and $p{$mm-1}[4] eq "+"){
				$seq=reverse substr($ccsseq,$p{$mm-1}[3],$k-$p{$mm-1}[3]-1);
				$seq=~tr/ATCGatcg/TAGCtacg/;
				$qv=reverse substr($ccsqv,$p{$mm-1}[3],$k-$p{$mm-1}[3]-1);
			}else{
				$seq=substr($ccsseq,$p{$mm-1}[3],$k-$p{$mm-1}[3]-1);
				$qv=substr($ccsqv,$p{$mm-1}[3],$k-$p{$mm-1}[3]-1);
			}

			$pp1="$p{$mm-1}[0]$p{$mm-1}[4]";
			$pp2="$p{$mm}[0]$p{$mm}[4]";
			$tag="UNKNOW";
			$tag="FL" if($pp1 eq "primer_S+" and $pp2 eq "primer_F-");
			$tag="FL" if($pp1 eq "primer_F+" and $pp2 eq "primer_S-");
			$tag="NFL" if($pp1 eq "primer_F-" and $pp2 eq "primer_S-");
			$tag="NFL" if($pp1 eq "primer_S+" and $pp2 eq "primer_F+");
			$tag="NFL" if($pp1 eq "primer_S-" and $pp2 eq "primer_S-");
			$tag="NFL" if($pp1 eq "primer_S+" and $pp2 eq "primer_S+");

			$barcode=reverse substr($seq,-16);
			$barcode=~tr/ATCGatcg/TAGCtacg/;			
			$umi=reverse substr($seq,-28,12);
			$umi=~tr/ATCGatcg/TAGCtacg/;
			
			$bcumi=reverse substr($seq,-28);
			$bcumi=~tr/ATCGatcg/TAGCtacg/;
			$bcumiqv=reverse substr($qv,-28);
			
			$dis=$k-$p{$mm-1}[3];
			if($tag eq "FL" and $sign{$mm-1}==0 and $dis>0){
				$tseq=substr($seq,0,-28);
				$tqv=substr($qv,0,-28);
				$pa=denovo_poly($tseq);
				$polyA="NA";$polyA=substr($tseq,$pa,);
				$Alen=0;$Alen=$polyA=~tr/A/A/;
				if($Alen>=3){
					$flseq=substr($tseq,0,$pa);
					$flqv=substr($qv,0,$pa);
					$palen=length($polyA);
				}else{
					$flseq=$tseq;
					$flqv=$tqv;
					$polyA="NA";$palen=0;
				}
				$seqlen=length($flseq);
				print FL2 "\@$tname\/$p{$mm-1}[3]\_$k\_CCS\n$flseq\n+\n$flqv\n" if($seqlen>=$min_seqlen);
				print BARCODE "$tname\/$p{$mm-1}[3]\_$k\_CCS,FL,$seqlen,$pp1,$p{$mm-1}[1],$p{$mm-1}[2],$pp2,$p{$mm}[1],$p{$mm}[2],$barcode,$umi,$palen,$polyA\n" if($seqlen>=$min_seqlen);

				print FL1 "\@$tname\/$p{$mm-1}[3]\_$k\_CCS\n$bcumi\n+\n$bcumiqv\n" if($seqlen>=$min_seqlen);
				
				$sign{$mm}=1;
			}elsif($tag eq "NFL"){
				$tseq=substr($seq,0,-28);$seqlen=length($tseq);
				$pa=denovo_poly($tseq);
				$polyA="NA";$polyA=substr($tseq,$pa,);
				$Alen=0;$Alen=$polyA=~tr/A/A/;
				if($Alen>=3){
					$flseq=substr($tseq,0,$pa);
					$palen=length($polyA);
				}else{
					$flseq=$tseq;
					$polyA="NA";$palen=0;
				}
				$seqlen=length($flseq);
				print BARCODE "$tname\/$p{$mm-1}[3]\_$k\_CCS,NFL,$seqlen,$pp1,$p{$mm-1}[1],$p{$mm-1}[2],$pp2,$p{$mm}[1],$p{$mm}[2],$barcode,$umi,$palen,$polyA\n" if($seqlen>=$min_seqlen);
			}elsif($tag eq "UNKNOW"){
				$tseq=$seq;$seqlen=length($tseq);
				print BARCODE "$tname\/$p{$mm-1}[3]\_$k\_CCS,UNKNOW,$seqlen,$pp1,$p{$mm-1}[1],$p{$mm-1}[2],$pp2,$p{$mm}[1],$p{$mm}[2],NA,NA,NA,NA\n" if($seqlen>=$min_seqlen);
			}
		}

		$mm++;
	}
}
close IN;
close FL1;
close FL2;
close BARCODE;

sub denovo_poly{
	$revseq=reverse($_[0]);$seqlen=length($revseq);
	$mark=0;%x=();
	while($revseq=~/(A+)/ig){
		$len=length $&;
		$pos=pos($revseq)-$len+1;
		$first=$pos if($mark==0);
		$end=$pos+$len;
		
		$x{$mark}[0]=$pos;$x{$mark}[1]=$end;
		$tail=substr($revseq,0,$end-1);
		$A=$tail=~tr/A/A/;$tAratio=0;$tAratio=sprintf("%.2f",$A/($end-$first)*100) if($end-1>0);
		
		$chunklen=$end-$x{$mark-1}[1] if($mark>0);
		$chunklen=$len if($mark==0);
		$chunk=substr($revseq,$x{$mark-1}[1]-1,$chunklen) if($mark>0);
		$chunk=substr($revseq,$pos-1,$chunklen) if($mark==0);
		$A=$chunk=~tr/A/A/;$cAratio=0;$cAratio=sprintf("%.2f",$A/$chunklen*100) if($chunklen>0);
		
		if($tAratio<70 or $cAratio<50){
			last;
		}
		$mark++;
	}
	$polystart=$seqlen-$x{$mark-1}[1]+1 if($mark>0);
	$polystart=$seqlen if($mark==0);
	return($polystart);
}
