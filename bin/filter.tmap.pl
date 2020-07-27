#! /usr/bin/perl

open IN,"$ARGV[0]";
%tmap;<IN>;
while(<IN>){
	chomp;
	@a=();@a=split();
	@b=();@b=split(/\./,$a[3]);
	$tmap{$b[0]}[0]++;
	$tmap{$b[0]}[1]=$_;
}
close IN;

@key1=();@key1=keys %tmap;
foreach $k(@key1){
	print "$tmap{$k}[1]\n" if($tmap{$k}[0]==1);
}
