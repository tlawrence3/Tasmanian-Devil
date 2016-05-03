#!/usr/bin/perl
use Math::Round;
use List::Util 'shuffle';

#Create files for the definitions of rules
open (NEGATIVECONTROLRULES, ">geneCalls_151015_negative_control_".$ARGV[0].".csv");
open (RULEBREAKDOWN, ">Rules_Negative_Control_".$ARGV[0].".txt");
open (PRINTLINES, ">Negative_Control_Raw_Analysis.csv");


#Read in the file with published iMM904 genes
open iMM904_Network_Genes, "<iMM904_Network_Genes.txt" or die$!;

while (<iMM904_Network_Genes>) {
	chomp;
	my ($iMM904_Network_Gene, $val1) = split /\t/;
	$iMM904{$iMM904_Network_Gene} = $val1;
	push (@iMM904_array,$iMM904_Network_Gene);
}
close iMM904_Network_Genes;

#Randomize the order of the array
@iMM904_shuffled = shuffle(@iMM904_array);

#Reverse the order of the array
my @iMM904_shuffled_reverse = reverse @iMM904_shuffled;

#From the command line input of fraction of genes and from the total number of iMM904 genes in the data set, determine the length of array 
my $count = keys %iMM904;
my $Microarray_Percent = $ARGV[0] * $count;
my $Microarray_Percent_Rounded = int($Microarray_Percent);

#Reduce the size of the array to the percentage of genes desired
splice @iMM904_shuffled, $Microarray_Percent_Rounded;

#Map the arrays to hashes to be able to look up the genes as keys later
%iMM904_shuffled_on = map { $_ => 1 } @iMM904_shuffled;

#Repeat the process for the reersed array
splice @iMM904_shuffled_reverse, $Microarray_Percent_Rounded;

%iMM904_shuffled_off = map { $_ => 1 } @iMM904_shuffled_reverse;


#Create rules from iMM904 genes for Ethanol
for $Genes1 ( keys %iMM904 ) {
#If gene is highly expressed, it is assigned a 1. If gene is lowly expressed, it is assigned a -1. 
	if (exists $iMM904_shuffled_on{$Genes1}) {
		push(@Negative_Control_Top,$Genes1);
		print NEGATIVECONTROLRULES "$Genes1,1\n";
	}
	elsif (exists $iMM904_shuffled_off{$Genes1}) {
		push(@Negative_Control_Bottom,$Genes1);
		print NEGATIVECONTROLRULES "$Genes1,-1\n";
	}
	#All other genes are assigend a 0
	else {
		push(@undefined_Negative_Control,$Genes1);
		print NEGATIVECONTROLRULES "$Genes1,0\n";
	}
} 

print RULEBREAKDOWN ("Highly expressed for Negative Control:", @Negative_Control_Top. "\n");
print RULEBREAKDOWN ("Lowly expressed for Negative Control:", @Negative_Control_Bottom. "\n");
print RULEBREAKDOWN ("Otherwise undefined for Negative Control:", @undefined_Negative_Control. "\n");
print RULEBREAKDOWN ("\n");

close (NEGATIVECONTROLRULES);
close (RULEBREAKDOWN);
close (PRINTLINES);
