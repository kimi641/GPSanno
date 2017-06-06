#!/usr/bin/perl -w
use strict;
use Text::CSV;
#perl script to pick useful information in HGMD csv files
sub loadHGMDSNP()
{
	my $hgmpSNP = shift;
	my %hash;
	my $csv = Text::CSV->new({sep_char => ',' });
	open(SNP, '<', $hgmpSNP) or die "No hgmpSNP csv";
	while(<SNP>){
		chomp;
		next if /^Type/;
		if($csv->parse($_)){
			my @fields	= $csv->fields();
			my $dbsnp	= $fields[3];
			my $coord	= $fields[5];
			my $hgvs	= $fields[6];
			my $gene	= $fields[8];
			my $disease	= $fields[9];
			my $pmid	= $fields[23];
			$hgvs =~ s/\.[0-9]+//;
			$hash{$hgvs} = [$gene, $disease, $pmid];
			print "$coord\t$dbsnp\t$hgvs\t$gene\t$disease\t$pmid\n";
		}else{
			#warn "Line could not be parsed:$_\n";
		}
	}
	close(SNP);
	return %hash;
	
}

sub loadHGMDindel()
{
	my $hgmdindel = shift;
	my %hash;
	my $csv = Text::CSV->new({sep_char => ',' });
	open(INDEL, '<', $hgmdindel) or die "No hgmpindel csv";
	while(<INDEL>){
		chomp;
		next if /^type/;
		if($csv->parse($_)){
			my @fields	= $csv->fields();
			my $gene	= $fields[2];
			my $disease	= $fields[3];
			my $hgvs	= $fields[4];
			my $coord	= $fields[6];
			my $pmid	= $fields[18];
			my $dbsnp	= $fields[20];
			$hgvs =~ s/\.[0-9]+//;
			$hash{$hgvs} = [$gene, $disease, $pmid];
			print "$coord\t$dbsnp\t$hgvs\t$gene\t$disease\t$pmid\n";
		}else{
			#warn "Line could not be parsed:$_\n";
		}
	}
	close(INDEL);
	return %hash;
}

sub main()
{
	print "#COORD\tDBSNP\tHGVS\tGENE\tDISEASE\tPMID\n";
	my %snphash = &loadHGMDSNP("HGMD_Advanced_Substitutions.csv");
	my %indelhash = &loadHGMDindel("HGMD_Advanced_Micro_Lesions.csv");
}

&main();
