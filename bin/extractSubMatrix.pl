use warnings;
use Getopt::Long;
use File::Basename;
use FindBin '$Bin';
use Cwd qw(abs_path);
sub usage {
    print STDERR<< "USAGE";

Despriprion: extract single cell submatrix of specified cell barcode list.
Author: shizhuoxing\@dingtalk.com
Data: 20210126
Usage: perl $0 -bclist bclist -inputematrix inputematrix -outputmatrix outputmatrix

Options:
        -bclist*:              specified cell barcode list
        -inputematrix*:        raw matrix in text format
        -outputmatrix*:        output submatrix in text format
        -help:                 print this help
USAGE
    exit 0;
}

GetOptions(
        "bclist:s" => \$bclist,
        "inputematrix:s" => \$inputematrix,
        "outputmatrix:s" => \$outputmatrix,
        "help:s" => \$help
);
die &usage() if ((!defined $bclist) or (!defined $inputematrix) or (!defined $outputmatrix) or (defined $help));

open IN,"$bclist";
%bc=();
while(<IN>){
	chomp;
	$bc{$_}=1;
}
close IN;

open IN,"$inputematrix";
open OUT,">$outputmatrix";
$head=<IN>;$head=~s/^\,//;
@a=();@a=split(/\,/,$head);
@pos=();$mark=0;$str="";
foreach $k(@a){
	$mark++;
	$k=~s/\"//g;$k=~s/\./\-/;
	push @pos,$mark if(defined $bc{$k});
	$str.="\"$k\"\," if(defined $bc{$k});
}

$str=~s/\,$//;
print OUT "$str\n";

while(<IN>){
	chomp;
	@a=();@a=split(/\,/,$_);
	print OUT "$a[0]";
	foreach $k(@pos){
		print OUT "\,$a[$k]";
	}
	print OUT "\n";
}
close IN;
close OUT;
