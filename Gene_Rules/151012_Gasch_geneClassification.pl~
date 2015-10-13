#!/usr/bin/perl
use Math::Round;

#Create files for the definitions of rules
open (ETHANOLRULES, ">geneCalls_151012_ethanol_".$ARGV[0].".csv");
open (GLUCOSERULES, ">geneCalls_151012_glucose_".$ARGV[0].".csv");
open (RULEBREAKDOWN, ">Rules_Gasch_".$ARGV[0].".txt");
open (PRINTLINES, ">Gasch_Raw_Analysis.csv");


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

#Read in microarray data for Ethanol condition from Gasch
open Microarray_Expression_Ethanol, "<151012_Gasch_ethanol.txt" or die$!;
while (<Microarray_Expression_Ethanol>) {
	chomp;
	my ($Microarray_Expression_Gene, $Microarray_Ethanol) = split /\t/;
	if (exists $iMM904{$Microarray_Expression_Gene}) {
		$Microarray_Exp_Ethanol{$Microarray_Expression_Gene} = $Microarray_Ethanol;
	}	
}
close Microarray_Expression_Ethanol;

#Read in microarray data for Glucose condition from Gasch
open Microarray_Expression_Glucose, "<151012_Gasch_glucose.txt" or die$!;
while (<Microarray_Expression_Glucose>) {
	chomp;
	my ($Microarray_Expression_Gene, $Microarray_Glucose) = split /\t/;
	if (exists $iMM904{$Microarray_Expression_Gene}) {
		$Microarray_Exp_Glucose{$Microarray_Expression_Gene} = $Microarray_Glucose;
	}	
}
close Microarray_Expression_Glucose;


#Create subroutines to sort the log averages for Glucose and Ethanol conditions
sub Microarray_Exp_Ethanol_Top {
	$Microarray_Exp_Ethanol{$b} <=> $Microarray_Exp_Ethanol{$a};
}

sub Microarray_Exp_Ethanol_Bottom {
	$Microarray_Exp_Ethanol{$a} <=> $Microarray_Exp_Ethanol{$b};
}

sub Microarray_Exp_Glucose_Top {
	$Microarray_Exp_Glucose{$b} <=> $Microarray_Exp_Glucose{$a};
}

sub Microarray_Exp_Glucose_Bottom {
	$Microarray_Exp_Glucose{$a} <=> $Microarray_Exp_Glucose{$b};
}


#Sort the log averages for a condiiton, then create an array with the sorted genes by value 
my @Microarray_Expression_Gene_Ethanol_Top_Sorted;
for my $Microarray_Expression_Gene_Ethanol_Top (sort Microarray_Exp_Ethanol_Top (keys(%Microarray_Exp_Ethanol))) {
	push(@Microarray_Expression_Gene_Ethanol_Top_Sorted,$Microarray_Expression_Gene_Ethanol_Top);
}

#From the command line input of fraction of genes and from the total number of iMM904 genes in the data set, determine the length of array 
my $count = keys %Microarray_Exp_Ethanol;
my $Microarray_Percent = $ARGV[0] * $count;
my $Microarray_Percent_Rounded = int($Microarray_Percent);

#Reduce the size of the array to the percentage of genes desired
splice @Microarray_Expression_Gene_Ethanol_Top_Sorted, $Microarray_Percent_Rounded;

#Map the arrays to hashes to be able to look up the genes as keys later
%Microarray_Ethanol_Top = map { $_ => 1 } @Microarray_Expression_Gene_Ethanol_Top_Sorted;

#Repeat for bottom expression
my @Microarray_Expression_Gene_Ethanol_Bottom_Sorted;
for my $Microarray_Expression_Gene_Ethanol_Bottom (sort Microarray_Exp_Ethanol_Bottom (keys(%Microarray_Exp_Ethanol))) {
	push(@Microarray_Expression_Gene_Ethanol_Bottom_Sorted,$Microarray_Expression_Gene_Ethanol_Bottom);
}

splice @Microarray_Expression_Gene_Ethanol_Bottom_Sorted, $Microarray_Percent_Rounded;

%Microarray_Ethanol_Bottom = map { $_ => 1 } @Microarray_Expression_Gene_Ethanol_Bottom_Sorted;

#Determine what percentage of expression the gene is expressed in for the iMM904 genes
$count_Ethanol = 0;
for (@Microarray_Expression_Gene_Ethanol_Bottom_Sorted) {
	$count_Ethanol += 1;
	#say $_;
	$count_div = int($count_Ethanol/$count*100)+1;
	$Microarray_Expression_Gene_Ethanol_Percent{$_} = $count_div;
}

#Repeat all operations for the Glucose data set
my @Microarray_Expression_Gene_Glucose_Top_Sorted;
for my $Microarray_Expression_Gene_Glucose_Top (sort Microarray_Exp_Glucose_Top (keys(%Microarray_Exp_Glucose))) {	
	push(@Microarray_Expression_Gene_Glucose_Top_Sorted,$Microarray_Expression_Gene_Glucose_Top);
}

splice @Microarray_Expression_Gene_Glucose_Top_Sorted, $Microarray_Percent_Rounded;

%Microarray_Glucose_Top = map { $_ => 1 } @Microarray_Expression_Gene_Glucose_Top_Sorted;

my @Microarray_Expression_Gene_Glucose_Bottom_Sorted;
for my $Microarray_Expression_Gene_Glucose_Bottom (sort Microarray_Exp_Glucose_Bottom (keys(%Microarray_Exp_Glucose))) {
	push(@Microarray_Expression_Gene_Glucose_Bottom_Sorted,$Microarray_Expression_Gene_Glucose_Bottom); 
}

splice @Microarray_Expression_Gene_Glucose_Bottom_Sorted, $Microarray_Percent_Rounded;

%Microarray_Glucose_Bottom = map { $_ => 1 } @Microarray_Expression_Gene_Glucose_Bottom_Sorted;

$count_Glucose = 0;
for (@Microarray_Expression_Gene_Glucose_Bottom_Sorted) {
	$count_Glucose += 1;
	#say $_;
	$count_div = int($count_Glucose/$count*100);
	$Microarray_Expression_Gene_Glucose_Percent{$_} = $count_div;
}

#Create rules from iMM904 genes for Ethanol
for $Genes1 ( keys %iMM904 ) {
#If gene is highly expressed, it is assigned a 1. If gene is lowly expressed, it is assigned a -1. 
	if (exists $Microarray_Ethanol_Top{$Genes1}) {
		push(@Ethanol_Top,$Genes1);
		print ETHANOLRULES "$Genes1,1\n";
	}
	elsif (exists $Microarray_Ethanol_Bottom{$Genes1}) {
		push(@Ethanol_Bottom,$Genes1);
		print ETHANOLRULES "$Genes1,-1\n";
	}
	#All other genes are assigend a 0
	else {
		push(@undefined_Ethanol,$Genes1);
		print ETHANOLRULES "$Genes1,0\n";
	}
} 

#Create rules from iMM904 genes for Glucose
for $Genes1 ( keys %iMM904 ) {
#If gene is highly expressed, it is assigned a 1. If gene is lowly expressed, it is assigned a -1. 
	if (exists $Microarray_Glucose_Top{$Genes1}) {
		push(@Glucose_Top,$Genes1);
		print GLUCOSERULES "$Genes1,1\n";
	}
	elsif (exists $Microarray_Glucose_Bottom{$Genes1}) {
		push(@Glucose_Bottom,$Genes1);
		print GLUCOSERULES "$Genes1,-1\n";
	}
	#All other genes are assigend a 0
	else {
		push(@undefined_Glucose,$Genes1);
		print GLUCOSERULES "$Genes1,0\n";
	}
} 

print RULEBREAKDOWN ("Highly expressed for Ethanol:", @Ethanol_Top. "\n");
print RULEBREAKDOWN ("Lowly expressed for Ethanol:", @Ethanol_Bottom. "\n");
print RULEBREAKDOWN ("Otherwise undefined Ethanol:", @undefined_Ethanol. "\n");
print RULEBREAKDOWN ("\n");

print RULEBREAKDOWN ("Highly expressed for Glucose:", @Glucose_Top. "\n");
print RULEBREAKDOWN ("Lowly expressed for Glucose:", @Glucose_Bottom. "\n");
print RULEBREAKDOWN ("Otherwise undefined Glucose:", @undefined_Glucose. "\n");
print RULEBREAKDOWN ("\n");

for $Genes1 ( keys %iMM904 ) { 
	print PRINTLINES ("$Genes1");
	print PRINTLINES (",");
	print PRINTLINES ("$Genes1");
	print PRINTLINES ("_");
	print PRINTLINES ("$Microarray_Exp_Ethanol{$Genes1}");
	print PRINTLINES ("_");
	print PRINTLINES ("$Microarray_Expression_Gene_Ethanol_Percent{$Genes1}");
	print PRINTLINES (",");
	print PRINTLINES ("$Genes1");
	print PRINTLINES ("_");
	print PRINTLINES ("$Microarray_Exp_Glucose{$Genes1}");
	print PRINTLINES ("_");
	print PRINTLINES ("$Microarray_Expression_Gene_Glucose_Percent{$Genes1}");
}


close (ETHANOLRULES);
close (GLUCOSERULES);
close (RULEBREAKDOWN);
close (PRINTLINES);
