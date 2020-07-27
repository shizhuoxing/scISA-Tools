#!/usr/bin/perl
use strict;
use Getopt::Std;
use FindBin;use lib $FindBin::Bin;

my $usage = q/Usage:
	sam2gff.pl [-R] [-o <output.gff3>] <input.sam> 
	Outputs a quick'n'dirty GFF3 conversion of spliced alignments given at stdin or as an input file name (SAM format expected) Option -R can be used to change the output format from GFF3 to the 	aligned segments as a simple "multi-region" format which is suitable for seqmanip/;

umask 0002;
getopts('Ro:') || die($usage."\n");
my %ids; #to avoid duplicates, need to keep track of all IDs $qname => HI
my $outfile=$Getopt::Std::opt_o;
if ($outfile) {
	open(OUTF, '>'.$outfile) || die("Error creating output file $outfile\n");
	select(OUTF);
}

# --
my $reg_out=$Getopt::Std::opt_R;
pop(@ARGV) while ($ARGV[-1] eq '-' || $ARGV[-1] eq 'stdin'); 
while (<>) {
	my $samline=$_;
	chomp;
	my ($qname, $flags, $gseq, $pos, $mapq, $cigar, $rnext, $pnext,$tlen, $seq, $qual, @extra)= split(/\t/);
	next if length($cigar)<2 || $pos==0;
	my $alnstrand= (($flags & 0x10)==0) ? '+' : '-';
	my @cigdata=($cigar=~m/(\d+[A-Z,=])/g);
	my $sflag = $flags & 0xc0;
	if ($sflag == 0x40) {
		$qname.='/1';
	} elsif ($sflag == 0x80) {
		$qname.='/2';
	}
	my ($ts, $xs);
	foreach my $tag (@extra) {
		if ($tag=~m/XS:A:([\+\-])/) {
			$xs=$1;
			last;
		}
		if ($tag=~m/ts:A:([\+\-])/) {
			$ts=$1;
			last;
		}
	}
	if ($ts && !$xs) {
		$xs=$ts;
		if (($flags & 0x10)!=0) {
			$xs = ($ts eq '-') ? '+' : '-';
		}
	}
 	my $tstrand=$xs || $alnstrand;
 	my ($mstart, $mend);
 	my @exons; #list of [exon_start, exon_end]
 	my $curpos=$pos;
 	$mstart=$pos;
 	foreach my $cd (@cigdata) {
		my $code=chop($cd);
 		#now $cd has the length of the CIGAR operation
		#next if $code eq 'S';
		if ($code eq 'M' || $code eq 'D') {
		#only advance genomic position for match and delete operations
			$curpos+=$cd;
 			$mend=$curpos-1; #advance the end of the exon
			next;
		}
		#if ($code eq 'I') {
		#   # skip reference gaps (extra read bases)
		#   next;
		#   $mend=$curpos-1;
		 #}
		if ($code eq 'N') { # intron
			#exon should really end here
			if ($mend) {
			push(@exons, [$mstart, $mend]);
			$mend=0;
			} else {
				die("Error: intron gap found before any exon?!\n$samline\n")
				unless @exons>0;
				# -- it can follow an Insertion, which followed another intron (!)
				# * in that case, the intron should be silendly extended, 
				#   and mstart (the next exon start) will be pushed accordingly
			}
			$curpos+=$cd;    # genomic position advancing to
			$mstart=$curpos; #    the start of the next exon
			$mend=0;
		}
	} #foreach cigar event
	#check the last interval
	if ($mend) {
		push(@exons, [$mstart, $mend]);
 	}
	my $hi=(++$ids{$qname}); #hit index
	$qname.=".m$hi";
	if ($reg_out) {
		print $qname."\t".$gseq.':'.join(',', (map { $$_[0].'-'.$$_[1] } @exons) );
		print '-' if $tstrand eq '-';
		print "\n";
		next;
	}
 	# print GFF3 here, as plain mRNA because it's easier
	print join("\t", $gseq, 'sam2gff', 'mRNA', $exons[0]->[0], 
	$exons[-1]->[1], '.', $tstrand, '.', "ID=$qname");
	#print additional tags as attributes?
	print "\n";
 	foreach my $exon (@exons) {
		print join("\t", $gseq, 'sam2gff', 'exon', $exon->[0], $exon->[1], '.', $tstrand, '.', "Parent=$qname\n");
	}
} #for each SAM line

# --
if ($outfile) {
	select(STDOUT);
	close(OUTF);
}

#************ Subroutines **************

sub checkOverlap { #check if segment $a-$b overlaps any of the exons in $rx
	my ($a, $b, $rx)=@_;
	return 0 if ($a>$$rx[-1]->[1] || $b<$$rx[0]->[0]); # not overlapping the transcript region at all
	foreach my $x (@$rx) {
		#check for exon overlap
		return 1 if ($a<=$$x[1] && $b>=$$x[0]);
		return 0 if $b<$$x[0]; #@$rx is sorted, so we can give up
	}
}
