#!/usr/bin/perl
use Math::Round;

#Create files for the definitions of rules
open (GCR1ARULES, ">geneCalls_141105_GCR1A_".$ARGV[0].".csv");
open (GCR1DRULES, ">geneCalls_141105_GCR1D_".$ARGV[0].".csv");
open (RULEBREAKDOWN, ">Rules_Baker_".$ARGV[0].".txt");
open (PRINTLINES, ">Baker_Raw_Analysis.csv");

#Read in files with differentially expressed genes for each condition
open Microarray_GCR1A_DE, '<141105_Baker_GCR1A_DE_genes.txt' or die$!;

while (<Microarray_GCR1A_DE>) {
	chomp;
	my ($Microarray_GCR1A_DE_Gene, $holder) = split /\t/;
	$GCR1A_DE_Genes{$Microarray_GCR1A_DE_Gene} = $holder;
}
close Microarray_GCR1A_DE;

open Microarray_GCR1D_DE, '<141105_Baker_GCR1D_DE_genes.txt' or die$!;

while (<Microarray_GCR1D_DE>) {
	chomp;
	my ($Microarray_GCR1D_DE_Gene, $holder) = split /\t/;
	$GCR1D_DE_Genes{$Microarray_GCR1D_DE_Gene} = $holder;
}
close Microarray_GCR1D_DE;

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

#Read in microarray data for GCR1 from Henry Baker with 3 samples per condition
open Microarray_Expression, "<Microarray_Expression.txt" or die$!;
while (<Microarray_Expression>) {
	chomp;
	my ($Microarray_Expression_Gene, $Microarray_GCR1A_1, $Microarray_GCR1A_2, $Microarray_GCR1A_3, $Microarray_GCR1D_1, $Microarray_GCR1D_2, $Microarray_GCR1D_3) = split /\t/;
	$Microarray_GCR1A_pre_log_avg{$Microarray_Expression_Gene} = int(($Microarray_GCR1A_1 + $Microarray_GCR1A_2 + $Microarray_GCR1A_3)/3); 
	$Microarray_GCR1D_pre_log_avg{$Microarray_Expression_Gene} = int(($Microarray_GCR1D_1 + $Microarray_GCR1D_2 + $Microarray_GCR1D_3)/3);
#Take log2 averages for each condition
	my $Microarray_GCR1A_1_log2 = logbase2($Microarray_GCR1A_1);
	my $Microarray_GCR1A_2_log2 = logbase2($Microarray_GCR1A_2);
	my $Microarray_GCR1A_3_log2 = logbase2($Microarray_GCR1A_3);
	my $Microarray_GCR1D_1_log2 = logbase2($Microarray_GCR1D_1);
	my $Microarray_GCR1D_2_log2 = logbase2($Microarray_GCR1D_2);
	my $Microarray_GCR1D_3_log2 = logbase2($Microarray_GCR1D_3);
	my $Microarray_GCR1A_avg = ($Microarray_GCR1A_1_log2 + $Microarray_GCR1A_2_log2 + $Microarray_GCR1A_3_log2)/2;
	my $Microarray_GCR1D_avg = ($Microarray_GCR1D_1_log2 + $Microarray_GCR1D_2_log2 + $Microarray_GCR1D_3_log2)/2;
	my $difflogavgsGCR1A = $Microarray_GCR1A_avg - $Microarray_GCR1D_avg;	
	my $difflogavgsGCR1D = $Microarray_GCR1D_avg - $Microarray_GCR1A_avg;
	$difflogavgsGCR1A_Microarray{$Microarray_Expression_Gene} = $difflogavgsGCR1A;
	$difflogavgsGCR1D_Microarray{$Microarray_Expression_Gene} = $difflogavgsGCR1D;
	$Microarray_Exp_GCR1A_1{$Microarray_Expression_Gene} = $Microarray_GCR1A_1;
	$Microarray_Exp_GCR1A_2{$Microarray_Expression_Gene} = $Microarray_GCR1A_2;
	$Microarray_Exp_GCR1A_3{$Microarray_Expression_Gene} = $Microarray_GCR1A_3;
	$Microarray_Exp_GCR1D_1{$Microarray_Expression_Gene} = $Microarray_GCR1D_1;
	$Microarray_Exp_GCR1D_2{$Microarray_Expression_Gene} = $Microarray_GCR1D_2;
	$Microarray_Exp_GCR1D_3{$Microarray_Expression_Gene} = $Microarray_GCR1D_3;
	$Microarray_Exp_GCR1A{$Microarray_Expression_Gene} = $Microarray_GCR1A_avg;
	$Microarray_Exp_GCR1D{$Microarray_Expression_Gene} = $Microarray_GCR1D_avg;	
}
close Microarray_Expression;

#Create subroutines to sort the log averages for GCR1A and GCR1D
sub Microarray_Exp_GCR1A_Top {
	$Microarray_Exp_GCR1A{$b} <=> $Microarray_Exp_GCR1A{$a};
}

sub Microarray_Exp_GCR1A_Bottom {
	$Microarray_Exp_GCR1A{$a} <=> $Microarray_Exp_GCR1A{$b};
}

sub Microarray_Exp_GCR1D_Top {
	$Microarray_Exp_GCR1D{$b} <=> $Microarray_Exp_GCR1D{$a};
}

sub Microarray_Exp_GCR1D_Bottom {
	$Microarray_Exp_GCR1D{$a} <=> $Microarray_Exp_GCR1D{$b};
}


#Sort the log averages for a condiiton, then create an array with the sorted genes by value 
my @Microarray_Expression_Gene_GCR1A_Top_Sorted;
for my $Microarray_Expression_Gene_GCR1A_Top (sort Microarray_Exp_GCR1A_Top (keys(%Microarray_Exp_GCR1A))) {
	push(@Microarray_Expression_Gene_GCR1A_Top_Sorted,$Microarray_Expression_Gene_GCR1A_Top);
}

#print "$Microarray_Expression_Gene_GCR1A_Top_Sorted\n";

#From the command line input of fraction of genes and from the total number of iMM904 genes in the data set, determine the length of array 
my $count = keys %Microarray_Exp_GCR1A;
my $Microarray_Percent = $ARGV[0] * $count;
my $Microarray_Percent_Rounded = int($Microarray_Percent);

#Reduce the size of the array to the percentage of genes desired
splice @Microarray_Expression_Gene_GCR1A_Top_Sorted, $Microarray_Percent_Rounded;

#Map the arrays to hashes to be able to look up the genes as keys later
%Microarray_GCR1A_Top = map { $_ => 1 } @Microarray_Expression_Gene_GCR1A_Top_Sorted;

my @Microarray_Expression_Gene_GCR1A_Bottom_Sorted;
for my $Microarray_Expression_Gene_GCR1A_Bottom (sort Microarray_Exp_GCR1A_Bottom (keys(%Microarray_Exp_GCR1A))) {
	push(@Microarray_Expression_Gene_GCR1A_Bottom_Sorted,$Microarray_Expression_Gene_GCR1A_Bottom);
}

#print $Microarray_Expression_Gene_GCR1A_Bottom_Sorted;

$count_GCR1A = 0;
for (@Microarray_Expression_Gene_GCR1A_Bottom_Sorted) {
	$count_GCR1A += 1;
	#say $_;
	$count_div = int($count_GCR1A/$count*100)+1;
	$Microarray_Expression_Gene_GCR1A_Percent{$_} = $count_div;
}

my @Microarray_Expression_Gene_GCR1D_Top_Sorted;
for my $Microarray_Expression_Gene_GCR1D_Top (sort Microarray_Exp_GCR1D_Top (keys(%Microarray_Exp_GCR1D))) {	
	push(@Microarray_Expression_Gene_GCR1D_Top_Sorted,$Microarray_Expression_Gene_GCR1D_Top);
}

splice @Microarray_Expression_Gene_GCR1D_Top_Sorted, $Microarray_Percent_Rounded;

%Microarray_GCR1D_Top = map { $_ => 1 } @Microarray_Expression_Gene_GCR1D_Top_Sorted;

my @Microarray_Expression_Gene_GCR1D_Bottom_Sorted;
for my $Microarray_Expression_Gene_GCR1D_Bottom (sort Microarray_Exp_GCR1D_Bottom (keys(%Microarray_Exp_GCR1D))) {
	push(@Microarray_Expression_Gene_GCR1D_Bottom_Sorted,$Microarray_Expression_Gene_GCR1D_Bottom); 
}

$count_GCR1D = 0;
for (@Microarray_Expression_Gene_GCR1D_Bottom_Sorted) {
	$count_GCR1D += 1;
	#say $_;
	$count_div = int($count_GCR1D/$count*100);
	$Microarray_Expression_Gene_GCR1D_Percent{$_} = $count_div;
}

#Create rules from iMM904 genes for GCR1A
for $Genes1 ( keys %iMM904 ) {
	#If gene is in differentially expressed genes for condition, it is assigned a 1
	if (exists $GCR1A_DE_Genes{$Genes1}) {
		push(@GCR1A_DE,$Genes1);
		print GCR1ARULES "$Genes1,1\n";
	}
#If gene is in differentially expressed genes for other condition, it is assigned a -1 unless it is in the top percentage of genes expressed, whereby it is assigned a 0
	elsif (exists $GCR1D_DE_Genes{$Genes1}) {
		if (exists $Microarray_GCR1A_Top{$Genes1}) {
			push(@GCR1A_DE_GCR1D_Top,$Genes1);
			print GCR1ARULES "$Genes1,0\n";
		}
		else {
			push(@GCR1A_DE_GCR1D,$Genes1);
			print GCR1ARULES "$Genes1,-1\n";
		}
	}
	#All other genes are assigend a 0
	else {
		push(@undefined_GCR1A,$Genes1);
		print GCR1ARULES "$Genes1,0\n";
	}
} 

#Repeat creating rules of GCR1D
for $Genes2 ( keys %iMM904 ) {
	#If gene is in differentially expressed genes for condition, it is assigned a 1
	if (exists $GCR1D_DE_Genes{$Genes2}) {
		push(@GCR1D_DE,$Genes2);
		print GCR1DRULES "$Genes2,1\n";
	}
#If gene is in differentially expressed genes for other condition, it is assigned a -1 unless it is in the top percentage of genes expressed, whereby it is assigned a 0
	elsif (exists $GCR1A_DE_Genes{$Genes2}) {
		if (exists $Microarray_GCR1D_Top{$Genes2}) {
			push(@GCR1D_DE_GCR1A_Top,$Genes2);
			print GCR1DRULES "$Genes2,0\n";
		}
		else {
			push(@GCR1D_DE_GCR1A,$Genes2);
			print GCR1DRULES "$Genes2,-1\n";
		}
	}
	#All other genes are assigend a 0
	else {
		push(@undefined_GCR1D,$Genes2);
		print GCR1DRULES "$Genes2,0\n";
	}
} 

print RULEBREAKDOWN ("Differentially expressed for GCR1A:", @GCR1A_DE. "\n");
print RULEBREAKDOWN ("Differentially expressed for GCR1D, in top percentage expressed for GCR1A:", @GCR1A_DE_GCR1D_Top. "\n");
print RULEBREAKDOWN ("Differentially expressed for GCR1D, not in top percentage expressed for GCR1A:", @GCR1A_DE_GCR1D."\n");	
print RULEBREAKDOWN ("Otherwise undefined GCR1A:", @undefined_GCR1A. "\n");
print RULEBREAKDOWN ("\n");

print RULEBREAKDOWN ("Differentially expressed for GCR1D:", @GCR1D_DE. "\n");
print RULEBREAKDOWN ("Differentially expressed for GCR1A, in top percentage expressed for GCR1D:", @GCR1D_DE_GCR1A_Top. "\n");
print RULEBREAKDOWN ("Differentially expressed for GCR1A, not in top percentage expressed for GCR1D:", @GCR1D_DE_GCR1A."\n");
print RULEBREAKDOWN ("Otherwise undefined GCR1D:", @undefined_GCR1D. "\n");		


for $Genes1 ( keys %iMM904 ) { 
	$difflogavgsGCR1A_Microarray_rounded = nearest(0.01,$difflogavgsGCR1A_Microarray{$Genes1});
	$difflogavgsGCR1D_Microarray_rounded = nearest(0.01,$difflogavgsGCR1D_Microarray{$Genes1});
#	print "$Genes1$Microarray_GCR1A_pre_log_avg{$Genes1},$Microarray_Expression_Gene_GCR1A_Percent{$Genes1},$difflogavgsGCR1A_Microarray_rounded,$Microarray_GCR1D_pre_log_avg{$Genes1},$Microarray_Expression_Gene_GCR1D_Percent{$Genes1},$difflogavgsGCR1D_Microarray_rounded\n";
	print PRINTLINES ("$Genes1");
	print PRINTLINES (",");
	print PRINTLINES ("$Genes1");
	print PRINTLINES ("_");
	print PRINTLINES ("$Microarray_GCR1A_pre_log_avg{$Genes1}");
	print PRINTLINES ("_");
	print PRINTLINES ("$Microarray_Expression_Gene_GCR1A_Percent{$Genes1}");
	print PRINTLINES ("_");
	print PRINTLINES ("$difflogavgsGCR1A_Microarray_rounded");
	print PRINTLINES (",");
	print PRINTLINES ("$Genes1");
	print PRINTLINES ("_");
	print PRINTLINES ("$Microarray_GCR1D_pre_log_avg{$Genes1}");
	print PRINTLINES ("_");
	print PRINTLINES ("$Microarray_Expression_Gene_GCR1D_Percent{$Genes1}");
	print PRINTLINES ("_");
	print PRINTLINES ("$difflogavgsGCR1D_Microarray_rounded\n");
}


close (GCR1ARULES);
close (GCR1DRULES);
close (RULEBREAKDOWN);
