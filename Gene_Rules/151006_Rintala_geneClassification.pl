#!/usr/bin/perl
use Math::Round;

#Create files for the definitions of rules
open (ANAEROBICRULES, ">geneCalls_151006_Anaerobic_".$ARGV[0].".csv");
open (AEROBICRULES, ">geneCalls_151006_Aerobic_".$ARGV[0].".csv");
open (RULEBREAKDOWN, ">Rules_Rintala_".$ARGV[0].".txt");
open (PRINTLINES, ">Rintala_Raw_Analysis.csv");


#Read in the file with published iMM904 genes
open iMM904_Network_Genes, "<iMM904_Network_Genes.txt" or die$!;

while (<iMM904_Network_Genes>) {
	chomp;
	my ($iMM904_Network_Gene, $val1) = split /\t/;
	$iMM904{$iMM904_Network_Gene} = $val1;
}
close iMM904_Network_Genes;

sub logbase2 {
	my $n = shift;
	return (log($n)/ log(2));
}

#Read in microarray data for anaerobic condition from Rintala
open Microarray_Expression_Anaerobic, "<151006_Rintala_0.txt" or die$!;
while (<Microarray_Expression_Anaerobic>) {
	chomp;
	my ($Microarray_Expression_Gene, $Microarray_Anaerobic) = split /\t/;
	$Microarray_Exp_Anaerobic{$Microarray_Expression_Gene} = $Microarray_Anaerobic;	
}
close Microarray_Expression_Anaerobic;

#Read in microarray data for aerobic condition from Rintala
open Microarray_Expression_Aerobic, "<151006_Rintala_20.9.txt" or die$!;
while (<Microarray_Expression_Aerobic>) {
	chomp;
	my ($Microarray_Expression_Gene, $Microarray_Aerobic) = split /\t/;
	$Microarray_Exp_Aerobic{$Microarray_Expression_Gene} = $Microarray_Aerobic;	
}
close Microarray_Expression_Aerobic;


#Create subroutines to sort the log averages for aerobic and anaerobic conditions
sub Microarray_Exp_Anaerobic_Top {
	$Microarray_Exp_Anaerobic{$b} <=> $Microarray_Exp_Anaerobic{$a};
}

sub Microarray_Exp_Anaerobic_Bottom {
	$Microarray_Exp_Anaerobic{$a} <=> $Microarray_Exp_Anaerobic{$b};
}

sub Microarray_Exp_Aerobic_Top {
	$Microarray_Exp_Aerobic{$b} <=> $Microarray_Exp_Aerobic{$a};
}

sub Microarray_Exp_Aerobic_Bottom {
	$Microarray_Exp_Aerobic{$a} <=> $Microarray_Exp_Aerobic{$b};
}


#Sort the log averages for a condiiton, then create an array with the sorted genes by value 
my @Microarray_Expression_Gene_Anaerobic_Top_Sorted;
for my $Microarray_Expression_Gene_Anaerobic_Top (sort Microarray_Exp_Anaerobic_Top (keys(%Microarray_Exp_Anaerobic))) {
	push(@Microarray_Expression_Gene_Anaerobic_Top_Sorted,$Microarray_Expression_Gene_Anaerobic_Top);
}

#From the command line input of fraction of genes and from the total number of iMM904 genes in the data set, determine the length of array 
my $count = keys %Microarray_Exp_Anaerobic;
my $Microarray_Percent = $ARGV[0] * $count;
my $Microarray_Percent_Rounded = int($Microarray_Percent);

#Reduce the size of the array to the percentage of genes desired
splice @Microarray_Expression_Gene_Anaerobic_Top_Sorted, $Microarray_Percent_Rounded;

#Map the arrays to hashes to be able to look up the genes as keys later
%Microarray_Anaerobic_Top = map { $_ => 1 } @Microarray_Expression_Gene_Anaerobic_Top_Sorted;

#Repeat for bottom expression
my @Microarray_Expression_Gene_Anaerobic_Bottom_Sorted;
for my $Microarray_Expression_Gene_Anaerobic_Bottom (sort Microarray_Exp_Anaerobic_Bottom (keys(%Microarray_Exp_Anaerobic))) {
	push(@Microarray_Expression_Gene_Anaerobic_Bottom_Sorted,$Microarray_Expression_Gene_Anaerobic_Bottom);
}

splice @Microarray_Expression_Gene_Anaerobic_Bottom_Sorted, $Microarray_Percent_Rounded;

%Microarray_Anaerobic_Bottom = map { $_ => 1 } @Microarray_Expression_Gene_Anaerobic_Bottom_Sorted;

#Determine what percentage of expression the gene is expressed in for the iMM904 genes
$count_Anaerobic = 0;
for (@Microarray_Expression_Gene_Anaerobic_Bottom_Sorted) {
	$count_Anaerobic += 1;
	#say $_;
	$count_div = int($count_Anaerobic/$count*100)+1;
	$Microarray_Expression_Gene_Anaerobic_Percent{$_} = $count_div;
}

#Repeat all operations for the aerobic data set
my @Microarray_Expression_Gene_Aerobic_Top_Sorted;
for my $Microarray_Expression_Gene_Aerobic_Top (sort Microarray_Exp_Aerobic_Top (keys(%Microarray_Exp_Aerobic))) {	
	push(@Microarray_Expression_Gene_Aerobic_Top_Sorted,$Microarray_Expression_Gene_Aerobic_Top);
}

splice @Microarray_Expression_Gene_Aerobic_Top_Sorted, $Microarray_Percent_Rounded;

%Microarray_Aerobic_Top = map { $_ => 1 } @Microarray_Expression_Gene_Aerobic_Top_Sorted;

my @Microarray_Expression_Gene_Aerobic_Bottom_Sorted;
for my $Microarray_Expression_Gene_Aerobic_Bottom (sort Microarray_Exp_Aerobic_Bottom (keys(%Microarray_Exp_Aerobic))) {
	push(@Microarray_Expression_Gene_Aerobic_Bottom_Sorted,$Microarray_Expression_Gene_Aerobic_Bottom); 
}

splice @Microarray_Expression_Gene_Aerobic_Bottom_Sorted, $Microarray_Percent_Rounded;

%Microarray_Aerobic_Bottom = map { $_ => 1 } @Microarray_Expression_Gene_Aerobic_Bottom_Sorted;

$count_Aerobic = 0;
for (@Microarray_Expression_Gene_Aerobic_Bottom_Sorted) {
	$count_Aerobic += 1;
	#say $_;
	$count_div = int($count_Aerobic/$count*100);
	$Microarray_Expression_Gene_Aerobic_Percent{$_} = $count_div;
}

#Create rules from iMM904 genes for Anaerobic
for $Genes1 ( keys %iMM904 ) {
#If gene is highly expressed, it is assigned a 1. If gene is lowly expressed, it is assigned a -1. 
	if (exists $Microarray_Anaerobic_Top{$Genes1}) {
		push(@Anaerobic_Top,$Genes1);
		print ANAEROBICRULES "$Genes1,1\n";
	}
	elsif (exists $Microarray_Anaerobic_Bottom{$Genes1}) {
		push(@Anaerobic_Bottom,$Genes1);
		print ANAEROBICRULES "$Genes1,-1\n";
	}
	#All other genes are assigend a 0
	else {
		push(@undefined_Anaerobic,$Genes1);
		print ANAEROBICRULES "$Genes1,0\n";
	}
} 

#Create rules from iMM904 genes for Aerobic
for $Genes1 ( keys %iMM904 ) {
#If gene is highly expressed, it is assigned a 1. If gene is lowly expressed, it is assigned a -1. 
	if (exists $Microarray_Aerobic_Top{$Genes1}) {
		push(@Aerobic_Top,$Genes1);
		print AEROBICRULES "$Genes1,1\n";
	}
	elsif (exists $Microarray_Aerobic_Bottom{$Genes1}) {
		push(@Aerobic_Bottom,$Genes1);
		print AEROBICRULES "$Genes1,-1\n";
	}
	#All other genes are assigend a 0
	else {
		push(@undefined_Aerobic,$Genes1);
		print AEROBICRULES "$Genes1,0\n";
	}
} 

print RULEBREAKDOWN ("Highly expressed for Anaerobic:", @Anaerobic_Top. "\n");
print RULEBREAKDOWN ("Lowly expressed for Anaerobic:", @Anaerobic_Bottom. "\n");
print RULEBREAKDOWN ("Otherwise undefined Anaerobic:", @undefined_Anaerobic. "\n");
print RULEBREAKDOWN ("\n");

print RULEBREAKDOWN ("Highly expressed for Aerobic:", @Aerobic_Top. "\n");
print RULEBREAKDOWN ("Lowly expressed for Aerobic:", @Aerobic_Bottom. "\n");
print RULEBREAKDOWN ("Otherwise undefined Aerobic:", @undefined_Aerobic. "\n");
print RULEBREAKDOWN ("\n");

for $Genes1 ( keys %iMM904 ) { 
	print PRINTLINES ("$Genes1");
	print PRINTLINES (",");
	print PRINTLINES ("$Genes1");
	print PRINTLINES ("_");
	print PRINTLINES ("$Microarray_Exp_Anaerobic{$Genes1}");
	print PRINTLINES ("_");
	print PRINTLINES ("$Microarray_Expression_Gene_Anaerobic_Percent{$Genes1}");
	print PRINTLINES (",");
	print PRINTLINES ("$Genes1");
	print PRINTLINES ("_");
	print PRINTLINES ("$Microarray_Exp_Aerobic{$Genes1}");
	print PRINTLINES ("_");
	print PRINTLINES ("$Microarray_Expression_Gene_Aerobic_Percent{$Genes1}");
}


close (ANAEROBICRULES);
close (AEROBICRULES);
close (RULEBREAKDOWN);
close (PRINTLINES);
