#!/usr/bin/perl -w
use strict;
use FindBin qw($Bin);
unshift (@INC, $Bin);
use LWP::Simple;
use utf8::all;
use JSON;
require "html.pm";
binmode(STDIN, ':encoding(utf8)');
binmode(STDOUT, ':encoding(utf8)');
binmode(STDERR, ':encoding(utf8)');

$| = 1;

our $percent = 0.0;
our $OBOPATH = "$Bin/dist/";
our $MAF_FILE = "$Bin/data/MAF.txt";
our $dbNSFP = "$Bin/data/dbNSFP.txt";
our %prediction;
our $table = "$Bin/data/case.tab";
our @report = ();
our $casehash;
our %hpohash;
our %keywords;
our $disease_hpo = "$Bin/data/disease_hpo.txt";
our $gene_disease = "$Bin/data/gene_disease.txt";
our %htmlhash = (
	"data"		=> "",
	"cmd"		=> "",
	"variants_all"	=> 0,
	"variants_fl"	=> 0,
	"variants_anno"	=> 0,
	"var_chr1"	=> 0,
	"var_chr2"	=> 0,
	"var_chr3"	=> 0,
	"var_chr4"	=> 0,
	"var_chr5"	=> 0,
	"var_chr6"	=> 0,
	"var_chr7"	=> 0,
	"var_chr8"	=> 0,
	"var_chr9"	=> 0,
	"var_chr10"	=> 0,
	"var_chr11"	=> 0,
	"var_chr12"	=> 0,
	"var_chr13"	=> 0,
	"var_chr14"	=> 0,
	"var_chr15"	=> 0,
	"var_chr16"	=> 0,
	"var_chr17"	=> 0,
	"var_chr18"	=> 0,
	"var_chr19"	=> 0,
	"var_chr20"	=> 0,
	"var_chr21"	=> 0,
	"var_chr22"	=> 0,
	"var_chrX"	=> 0,
	"var_chrY"	=> 0,
	"var_chr1_an"	=> 0,
	"var_chr2_an"	=> 0,
	"var_chr3_an"	=> 0,
	"var_chr4_an"	=> 0,
	"var_chr5_an"	=> 0,
	"var_chr6_an"	=> 0,
	"var_chr7_an"	=> 0,
	"var_chr8_an"	=> 0,
	"var_chr9_an"	=> 0,
	"var_chr10_an"	=> 0,
	"var_chr11_an"	=> 0,
	"var_chr12_an"	=> 0,
	"var_chr13_an"	=> 0,
	"var_chr14_an"	=> 0,
	"var_chr15_an"	=> 0,
	"var_chr16_an"	=> 0,
	"var_chr17_an"	=> 0,
	"var_chr18_an"	=> 0,
	"var_chr19_an"	=> 0,
	"var_chr20_an"	=> 0,
	"var_chr21_an"	=> 0,
	"var_chr22_an"	=> 0,
	"var_chrX_an"	=> 0,
	"var_chrY_an"	=> 0,
	"frameshift_deletion"		=> 0,
	"frameshift_insertion"		=> 0,
	"nonframeshift_deletion"	=> 0,
	"nonframeshift_insertion"	=> 0,
	"nonsynonymous_SNV"		=> 0,
	"stopgain"			=> 0,
	"stoploss"			=> 0,
	"synonymous_SNV"		=> 0,
	"unknown"			=> 0,
	"probably_damaging"	=> 0,
	"possibly_damaging"	=> 0,
	"benign"		=> 0
	);

#Use HGMD annotate annovar results

my ($vcfin, $out, $maf, $family) = &getopts(@ARGV);
chop($out);
our $out_path = substr($out,0, rindex($out,"\/"));
chdir($out_path);
our $log = "gps.log";
our $annovarout = "temp";
our $avi_file = "temp.avinput";
our $ffin = "temp.family.vcf";
our $filter_vcf = "temp.filter.vcf";
our $OBOOUT = "tmp/annotationsCONCEPTS.txt";
open (LOG,'>>',$log) or die "Can't open $log:$!";

#Show help
sub help()
{
	print "Usage:\n";
	print "\tGPSanno.pl [options] -i input.vcf -o outputdir\n";
	printf "\t%-4s:\t%s\n","-i","input VCF file";
	printf "\t%-4s:\t%s\n","-o","output directory name";
	printf "\t%-4s:\t%s\n","-maf","MAF cutoff(default 0.01)";
	printf "\t%-4s:\t%s\n","-p","Patient tags in merge vcf eg:child,brother,sister";
	printf "\t%-4s:\t%s\n","-n","Normal tags in merge vcf eg:father,mother";
	print "Example:\n";
	print "\t./GPSanno.pl -i example.vcf -o example [-maf 0.01] [-p child] [-n father,mother]\n";
}

#get options
sub getopts()
{
	my @argv = @_;
	&help && exit(0) unless $#argv > 0;
	my $in;
	my $out;
	my $pph = "";
	my $maf = 0.01;
	my $patient = '';
	my $normal = '';
	while(my $argc = (shift @argv)){
		if($argc eq "-i"){
			$in = (shift @argv);
		}elsif($argc eq "-o"){
			$out = (shift @argv);
		}elsif($argc eq "-maf"){
			$maf = (shift @argv);
		}elsif($argc eq "-p"){
			$patient = (shift @argv);
		}elsif($argc eq "-n"){
			$normal = (shift @argv);
		}elsif($argc eq "-h"){
			&help;
			exit(0);
		}else{
			&help;
			exit(1);
		}
	}
	return $in, $out, $maf, $patient, $normal;
}

#get Abstract from pubmed
sub getAbstract()
{
	my $pubmedID = shift;
	my $url = "";
	$url .= get( "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=$pubmedID&rettype=medline&retmode=xml" );
	my $abstract = "";
	#warn "start get abstract";
	$abstract .= $1 while ($url =~ s/<AbstractText.*>(.*)<\/AbstractText.*>//);
	my $keywordslist = "";
	$keywordslist = $1 if ($url =~ /<KeywordList.*>([\s\S]*)<\/KeywordList>/);
	while ($keywordslist =~ /<Keyword.*>(.*)<\/Keyword>/){
		my $keyword = $1;
		if (exists $keywords{$keyword}){
			$keywords{$keyword} += 1;
		}else{
			$keywords{$keyword} = 1;
		}
		$keywordslist =~ s/<Keyword.*>(.*)<\/Keyword>//;
	}
	my $meshlist = "";
	$meshlist = $1 if ($url =~ /<MeshHeadingList>([\s\S]*)<\/MeshHeadingList>/);
	while ($meshlist =~ /<DescriptorName.*(.*)<\/<DescriptorName>/){
		my $mesh = $1;
		$abstract .= $mesh;
		if (exists $keywords{$mesh}){
			$keywords{$mesh} += 1;
		}else{
			$keywords{$mesh} = 1;
		}
		$meshlist =~ s/<DescriptorName.*(.*)<\/<DescriptorName>//;
	}
	#warn "abstract\n";
	return $abstract;
}

#remove duplicated variables from array
sub uniq{
	my %seen;
	return grep { !$seen{$_}++ }@_;
}

#Use OBO-annotator search HPO from astract
sub getOBO(){
	my $pmid = shift;
	my $abstract = shift;
	my $input = "input.txt";
	my $outdir = "tmp/";
	open (IN, '>', $input) or die "Can't open abstract.txt";
	print IN $pmid,"\n";
	print IN $abstract;
	close(IN);
	mkdir ("tmp");
	my $obojar = $OBOPATH . "OBOAnnotatorNoGUI.jar";
	`java -jar $obojar $input $outdir`;
	my @hpos;
	my $str = `cut -f 2 $OBOOUT | uniq`;
	my @ids = split(/\n/, $str);
	shift @ids;
	shift @ids;
	push @hpos, @ids;
	@hpos = &uniq(@hpos);
	unlink $input;
	`rm -r $outdir`;
	return @hpos
}

#Load case datasets
sub loadCASE()
{
	my $casefile = shift;
	my %hash;
	open(IN, '<', $casefile) or die "Can't open $casefile:$!";
	while(<IN>){
		chomp;
		next if (/^#/);
		my @fields = split(/\t/);
		my ($coord, $dbsnp, $hgvs, $gene, $disease, $pmid )= @fields;
		if (exists $hash{$gene}){
			if($hgvs ne "NULL"){
				if (!$hash{$gene}){
					my %new_hash = ($hgvs => [$disease, $pmid]);
					$hash{$gene} = \%new_hash;
				}else{
					$hash{$gene}{$hgvs} = [$disease, $pmid];
				}
			}
		}else{
			if($hgvs ne "NULL"){
				my %new_hash = ($hgvs => [$disease, $pmid]);
				$hash{$gene} = \%new_hash;
			}
		}
	}		
	close(IN);
	return \%hash;
}


#Get HPO from medgen databse
sub medgen()
{
	my $gene = shift;
	my $disease = shift;
	my $xml = get( "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=medgen&term=$disease+AND+$gene" . "[Gene]&retmax=10&sort=relevance&usehistory=y" );
	my @medgenid;
	my @hpos;
	return unless defined $xml;
	push( @medgenid, $1) while ($xml =~ s/<Id>(.*)<\/Id>//);
	for my $id (@medgenid){
		my $mxml = get( "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=medgen&id=$id&retmode=xml" );
		my $clinicalfs = $1 if ($mxml =~ /<ClinicalFeatures>(.*)<\/ClinicalFeatures>/);
		push(@hpos, $1) while ($clinicalfs =~ s/<Name>([\w|\s]+)<\/Name>//);
	}
	return &uniq(@hpos);
}

#Use e-tuils search pubmed including genes
sub pmidFromgene()
{
	my $gene = shift;
	my $disease = shift;
	my $hpo = shift;
	my @pmids;
	my $xml = get( "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=$gene" . "[Title/Abstract]+AND+$disease" . "[Title/Abstract]" . "+AND+$hpo&retmax=10&sort=relevance&usehistory=y" );
	if (defined $xml){
		push( @pmids, $1 ) while ($xml =~ s/<Id>(.*)<\/Id>//);
		return &uniq(@pmids);
	}else{
		return;
	}
}

sub searchpre
{
	my $coord = shift;
	my $alt = shift;
	$coord =~ s/chr// if $coord =~ /chr/;
	my $prediction = "Unknown";
	my $result = `grep $coord $dbNSFP`;
	chomp($result);
	if ($result){
		my @cols = split(/\t/,$result);
		my @alts = split(",",$cols[2]);
		my @scores = split(",",$cols[3]);
		for (my $i =0; $i <=$#alts;$i++){
			if ($alt = $alts[$i]){
				if($scores[$i] =~ /D/){
					$prediction = "probably_damaging";
				}elsif($scores[$i] =~ /P/){
    					$prediction = "possibly_damaging";
				}elsif($scores[$i] =~ /B/){
        				$prediction = "benign";
				}
			}
		}
		
	}
	return $prediction;
}

#search gene_disease
sub searchdiseases()
{
	my $gene = shift;
	my $lines = `grep -i "\\\\b$gene\\\\b" $gene_disease|cut -f 2`;
	my @diseases = split(/\n/, $lines);
	return @diseases;
}
sub casesearch()
{
	my $gene = shift;
	my $hgvs = shift;
	my $type = shift;
	my $prediction = shift;
	my $pmid = "null";
	my $count = 0;
	my $disease = "null";
	if (exists $$casehash{$gene}){
		my $genehash = $$casehash{$gene};
		if (exists $$genehash{$hgvs}){
			$disease = $$genehash{$hgvs}[0];
			$pmid = $$genehash{$hgvs}[1];
			if ($pmid =~ /[0-9]+/ ){
				my $abstract = &getAbstract($pmid);
				my @hpos = &getOBO($pmid, $abstract);
				foreach(@hpos){
					my $key = "$_\t$disease\t$gene\t$type\t$prediction";
					if (exists $hpohash{$key}){
						my $subhash = $hpohash{$key};
						$$subhash{"HGMD"} += 1;
						$count += 1;
					}else{
						my %subhash = (
							"HGMD"		=> 0,
							"GENE_MED"	=> 0,
							"GENE_HPO"	=> 0,
						);
						$subhash{"HGMD"} += 1;
						$count += 1;
						$hpohash{$key} = \%subhash;
					}
				}
			}
			my $disease_gene_records = `grep -i "$disease" $disease_hpo|cut -f 3`;
			my @hpos = split(/\n/, $disease_gene_records);
			foreach(@hpos){
				my $key = "$_\t$disease\t$gene\t$type\t$prediction";
				if (exists $hpohash{$key}){
					my $subhash = $hpohash{$key};
					$$subhash{"HGMD"} += 1;
					$count += 1;
				}else{
					my %subhash = (
						"HGMD"		=> 0,
						"GENE_MED"	=> 0,
						"GENE_HPO"	=> 0,
					);
					$subhash{"HGMD"} += 1;
					$count += 1;
					$hpohash{$key} = \%subhash;
				}
			}
		}
	}
	#warn "$pmcount\n";
	return ($count,$disease, $pmid);
}

sub medsearch{
	my $gene = shift;
	my $type = shift;
	my $prediction = shift;
	my $count = 0;
	my @diseases = &searchdiseases($gene);
	if (@diseases){
		for my $disease(@diseases){
			my @hpos = &medgen($gene, $disease);
			foreach(@hpos){
				my $key = "$_\t$disease\t$gene\t$type\t$prediction";
				if (exists $hpohash{$key}){
					my $subhash = $hpohash{$key};
					$$subhash{"GENE_MED"} += 1;
					$count += 1;
				}else{
					my %subhash = (
						"HGMD"		=> 0,
						"GENE_MED"	=> 0,
						"GENE_HPO"	=> 0,
					);
					$subhash{"GENE_MED"} += 1;
					$hpohash{$key} = \%subhash;
					$count += 1;
				}
			}	
		}
	}else{
		return 0;
	}
	return $count;
}
sub hposearch{
	my $gene = shift;
	my $type = shift;
	my $prediction = shift;
	my @pmids = ();
	my $count = 0;
	my @diseases = &searchdiseases($gene);
	if (@diseases){
		for my $disease (@diseases){
			my $records = `grep -i "$disease" $disease_hpo|cut -f 3`;
			my @hpos = split(/\n/, $records);
			@hpos = &uniq(@hpos);
			for my $hpo (@hpos){
				my $key = "$hpo\t$disease\t$gene\t$type\t$prediction";
				if (exists $hpohash{$key}){
					my $subhash = $hpohash{$key};
					$$subhash{"GENE_HPO"} += 1;
					$count += 1;
				}else{
					my %subhash = (
						"HGMD"		=> 0,
						"GENE_MED"	=> 0,
						"GENE_HPO"	=> 0,
					);
					$subhash{"GENE_HPO"} += 1;
					$hpohash{$key} = \%subhash;
					$count += 1;
				}
				my @thispm = &pmidFromgene($disease, $gene, $hpo);
				push (@pmids, @thispm);
			}	
		}
	}else{
		return 0;
	}
	@pmids = &uniq(@pmids);
	return $count, \@pmids;
}
#annotate ANNOVAR results
sub anno()
{
	my $annovarout = shift;
	if (-e "$annovarout.exonic_variant_function"){
		open (EVF, '<', "$annovarout.exonic_variant_function") or die "No annovar output\n";
	}
	my $now = 0;
	my $total = `wc -l "$annovarout.exonic_variant_function" | cut -d ' ' -f 1`;
	chomp($total);
	while(<EVF>)
	{
		chomp;
		my @fields	= split(/\t/);
		my $line	= $fields[0];
		my $type	= $fields[1];
		my $info_hgvs 	= $fields[2];
		my $chr		= $fields[3];
		my $start	= $fields[4];
		my $end		= $fields[5];
		my $ref		= $fields[6];
		my $alt		= $fields[7];
		my @hpocounts	= (0,0,0);	#hpocounts[HGMD_PMID,HGMD_HPO,MEDGEN,HPO]
		my @diseases;
		my $aaref	= "";
		my $aaalt	= "";
		my @pm;
		$now += 1;
		print "\r";
		print "$now/$total";
		$percent += (1/$total) * 0.8;
		printf LOG "%d%%\n",$percent * 100;
		$htmlhash{"variants_fl"} += 1;
		if($chr eq "chr1" or $chr eq "1"){
			$htmlhash{"var_chr1"} += 1;
		}elsif($chr eq "chr2" or $chr eq "2"){
			$htmlhash{"var_chr2"} += 1;
		}elsif($chr eq "chr3" or $chr eq "3"){
			$htmlhash{"var_chr3"} += 1;
		}elsif($chr eq "chr4" or $chr eq "4"){
			$htmlhash{"var_chr4"} += 1;
		}elsif($chr eq "chr5" or $chr eq "5"){
			$htmlhash{"var_chr5"} += 1;
		}elsif($chr eq "chr6" or $chr eq "6"){
			$htmlhash{"var_chr6"} += 1;
		}elsif($chr eq "chr7" or $chr eq "7"){
			$htmlhash{"var_chr7"} += 1;
		}elsif($chr eq "chr8" or $chr eq "8"){
			$htmlhash{"var_chr8"} += 1;
		}elsif($chr eq "chr9" or $chr eq "9"){
			$htmlhash{"var_chr9"} += 1;
		}elsif($chr eq "chr10" or $chr eq "10"){
			$htmlhash{"var_chr10"} += 1;
		}elsif($chr eq "chr11" or $chr eq "11"){
			$htmlhash{"var_chr11"} += 1;
		}elsif($chr eq "chr12" or $chr eq "12"){
			$htmlhash{"var_chr12"} += 1;
		}elsif($chr eq "chr13" or $chr eq "13"){
			$htmlhash{"var_chr13"} += 1;
		}elsif($chr eq "chr14" or $chr eq "14"){
			$htmlhash{"var_chr14"} += 1;
		}elsif($chr eq "chr15" or $chr eq "15"){
			$htmlhash{"var_chr15"} += 1;
		}elsif($chr eq "chr16" or $chr eq "16"){
			$htmlhash{"var_chr16"} += 1;
		}elsif($chr eq "chr17" or $chr eq "17"){
			$htmlhash{"var_chr17"} += 1;
		}elsif($chr eq "chr18" or $chr eq "18"){
			$htmlhash{"var_chr18"} += 1;
		}elsif($chr eq "chr19" or $chr eq "19"){
			$htmlhash{"var_chr19"} += 1;
		}elsif($chr eq "chr20" or $chr eq "21"){
			$htmlhash{"var_chr20"} += 1;
		}elsif($chr eq "chr21" or $chr eq "21"){
			$htmlhash{"var_chr21"} += 1;
		}elsif($chr eq "chr22" or $chr eq "22"){
			$htmlhash{"var_chr22"} += 1;
		}elsif($chr eq "chrX" or $chr eq "X"){
			$htmlhash{"var_chrX"} += 1;
		}elsif($chr eq "chrY" or $chr eq "Y"){
			$htmlhash{"var_chrY"} += 1;
		}
		my $typekey = $type;
		$typekey  =~ s/\s/_/g;
		$htmlhash{$typekey} += 1;
		next if $info_hgvs =~ /UNKNOWN/;
		next if $type eq "synonymous SNV"; 
		#get pph2 prediction
		my $pph2_key = "$chr\_$start";
		my $pph2_value;
		my $precoord = "$chr"."_$start";
		my $prediction = &searchpre($precoord,$alt);
		$htmlhash{$prediction} += 1;
		#next if $prediction =~ /benign/;
		if ($type =~ /frameshift/){
			$prediction = $type;
		}elsif($type =~ /stop/){
			$prediction = $type;
		}
		my @hgvss;
		my @genes;
		my @infolist = split(/,/,$info_hgvs);
		#Get gene and hgvs info
		for my $info (@infolist){
			my @eles = split(/:/,$info);
			my $gene = $eles[0];
			my $nmid = $eles[1];
			my $coord;
			my $hgvs;
			if (@eles > 3){
				my $nas = $eles[3];
				$nas =~ s/>/-/;
				if ($nas =~ /ins/){
					my $aas = $eles[4];
					if($aas =~ /p\.([A-Za-z]+)\d+delins([A-Za-z]+)/){
						$aaref = $1;
						$aaalt = $2;
					}
				}elsif($nas =~ /del/){
					$aaref = "";
					$aaalt = "";
				}else{
					my $aas = $eles[4];
					if ($aas =~ /p\.([A-Za-z]+)\d+([A-Za-z]+)/){
						$aaref = $1;
						$aaalt = $2;
					}
				}
				$hgvs = "$nmid:$nas";
				push(@hgvss, $hgvs);
			}else{
				$hgvs = $nmid;
				push(@hgvss, $hgvs);
			}
			push (@genes,$gene);
			my ($hgcount, $disease, $pmid) = &casesearch($gene,$hgvs,$type,$prediction);
			$hpocounts[0] += $hgcount;
			push (@diseases, $disease) unless $disease =~ /null/;
			push (@pm, $pmid) unless $pmid =~ /null/;
			#warn $hpocounts[0],"\n";
		}
		@genes = &uniq(@genes);
		if ($hpocounts[0] > 0){
			goto PTABLE;
		}
		for my $gene(@genes){
			my $medcount = &medsearch($gene,$type,$prediction);
			$hpocounts[1] += $medcount;
			if ($medcount != 0){
				push(@diseases, &searchdiseases($gene));
				next;
			}
			my ($hpocount, $ppm)  = &hposearch($gene,$type,$prediction);
			$hpocounts[2] += $hpocount;
			if (defined $ppm){
				push(@pm, @$ppm);
			}
			push(@diseases, &searchdiseases($gene));
		}
		PTABLE:
		if ($type =~ /deletion/){
			$start -= 1;
		}
		my $str = "$chr\t$start\t$ref\t$alt\t";
		@diseases = &uniq(@diseases);
		foreach(@diseases){
			$str .= "$_,";
		}
		chop($str) if $str =~ /,$/;
		$str .= "\t";
		foreach(@genes){
			$str .= "$_,";
		}
		chop($str) if $str =~ /,$/;
		$str .= "\t$type\t$aaref\t$aaalt\t$prediction\t";
		#warn $hpocounts[0],"\n";
		$str .= "$hpocounts[0]\t$hpocounts[1]\t$hpocounts[2]\t";
		@pm = &uniq(@pm);
		foreach(@pm){
			$str .= "$_,";
		}
		chop($str) if $str =~ /,$/;
		$str .= "\n";
		push (@report, $str);
		#warn "one loop\n";
	}
	close(EVF);
	if (-e "$annovarout.variant_function"){
		open (VF, '<', "$annovarout.variant_function") or die "No annovar output\n";
	}
	while(<VF>)
	{
		chomp;
		next unless /^splicing/;
		my @fields	= split(/\t/);
		my $type	= $fields[0];
		my $info_hgvs	= $fields[1];
		my $chr 	= $fields[2];
		my $start	= $fields[3];
		my $end		= $fields[4];
		my $ref		= $fields[5];
		my $alt		= $fields[6];
		my @hpocounts	= (0,0,0);	#hpocounts[HGMD_PMID,HGMD_HPO,MEDGEN,HPO]
		my @diseases;
		my $aaref	= "";
		my $aaalt	= "";
		my @pm;
		$htmlhash{"variants_fl"} += 1;
		if($chr eq "chr1" or $chr eq "1"){
			$htmlhash{"var_chr1"} += 1;
		}elsif($chr eq "chr2" or $chr eq "2"){
			$htmlhash{"var_chr2"} += 1;
		}elsif($chr eq "chr3" or $chr eq "3"){
			$htmlhash{"var_chr3"} += 1;
		}elsif($chr eq "chr4" or $chr eq "4"){
			$htmlhash{"var_chr4"} += 1;
		}elsif($chr eq "chr5" or $chr eq "5"){
			$htmlhash{"var_chr5"} += 1;
		}elsif($chr eq "chr6" or $chr eq "6"){
			$htmlhash{"var_chr6"} += 1;
		}elsif($chr eq "chr7" or $chr eq "7"){
			$htmlhash{"var_chr7"} += 1;
		}elsif($chr eq "chr8" or $chr eq "8"){
			$htmlhash{"var_chr8"} += 1;
		}elsif($chr eq "chr9" or $chr eq "9"){
			$htmlhash{"var_chr9"} += 1;
		}elsif($chr eq "chr10" or $chr eq "10"){
			$htmlhash{"var_chr10"} += 1;
		}elsif($chr eq "chr11" or $chr eq "11"){
			$htmlhash{"var_chr11"} += 1;
		}elsif($chr eq "chr12" or $chr eq "12"){
			$htmlhash{"var_chr12"} += 1;
		}elsif($chr eq "chr13" or $chr eq "13"){
			$htmlhash{"var_chr13"} += 1;
		}elsif($chr eq "chr14" or $chr eq "14"){
			$htmlhash{"var_chr14"} += 1;
		}elsif($chr eq "chr15" or $chr eq "15"){
			$htmlhash{"var_chr15"} += 1;
		}elsif($chr eq "chr16" or $chr eq "16"){
			$htmlhash{"var_chr16"} += 1;
		}elsif($chr eq "chr17" or $chr eq "17"){
			$htmlhash{"var_chr17"} += 1;
		}elsif($chr eq "chr18" or $chr eq "18"){
			$htmlhash{"var_chr18"} += 1;
		}elsif($chr eq "chr19" or $chr eq "19"){
			$htmlhash{"var_chr19"} += 1;
		}elsif($chr eq "chr20" or $chr eq "21"){
			$htmlhash{"var_chr20"} += 1;
		}elsif($chr eq "chr21" or $chr eq "21"){
			$htmlhash{"var_chr21"} += 1;
		}elsif($chr eq "chr22" or $chr eq "22"){
			$htmlhash{"var_chr22"} += 1;
		}elsif($chr eq "chrX" or $chr eq "X"){
			$htmlhash{"var_chrX"} += 1;
		}elsif($chr eq "chrY" or $chr eq "Y"){
			$htmlhash{"var_chrY"} += 1;
		}
		my $typekey = $type;
		$typekey  =~ s/\s/_/g;
		$htmlhash{$typekey} += 1;
		next if $info_hgvs =~ /UNKNOWN/;
		next if $type eq "synonymous SNV"; 
		#get pph2 prediction
		my $pph2_key = "$chr\_$start";
		my $pph2_value;
		my $precoord = "$chr"."_$start";
		my $prediction = &searchpre($precoord,$alt);
		$htmlhash{$prediction} += 1;
		#next if $prediction =~ /benign/;
		if ($type =~ /frameshift/){
			$prediction = $type;
		}elsif($type =~ /stop/){
			$prediction = $type;
		}
		my @hgvss;
		my $gene;
		my $hgvs_info;
		if ($info_hgvs =~ /^(\w+)$/){
			$gene = $1;
		}
		elsif ($info_hgvs =~ /^(\w+)\((.+)\)/){
			$gene = $1;
			$hgvs_info = $2;
		}
		if (defined $hgvs_info){
			my @infolist = split(/,/,$hgvs_info);
			#Get gene and hgvs info
			for my $info (@infolist){
				my @eles = split(/:/,$info);
				my $nmid = $eles[0];
				my $coord;
				my $hgvs;
				if (@eles > 2){
					my $nas = $eles[2];
					$nas =~ s/>/-/;
					$aaref = "";
					$aaalt = "";
					$hgvs = "$nmid:$nas";
					push(@hgvss, $hgvs);
				}else{
					$hgvs = $nmid;
					push(@hgvss, $hgvs);
				}
				my ($hgcount, $disease, $pmid) = &casesearch($gene,$hgvs,$type,$prediction);
				$hpocounts[0] += $hgcount;
				push (@diseases, $disease) unless $disease =~ /null/;
				push (@pm, $pmid) unless $pmid =~ /null/;
				#warn $hpocounts[0],"\n";
			}
		}
		if ($hpocounts[0] > 0){
			goto PTABLE;
		}
		my $medcount = &medsearch($gene,$type,$prediction);
		$hpocounts[1] += $medcount;
		if ($medcount != 0){
			push(@diseases, &searchdiseases($gene));
			goto PTABLE;
		}
		my ($hpocount, $ppm)  = &hposearch($gene,$type,$prediction);
		$hpocounts[2] += $hpocount;
		if (defined $ppm){
			push(@pm, @$ppm);
		}
		push(@diseases, &searchdiseases($gene));
		PTABLE:
		my $str = "$chr\t$start\t$ref\t$alt\t";
		@diseases = &uniq(@diseases);
		foreach(@diseases){
			$str .= "$_,";
		}
		chop($str) if $str =~ /,$/;
		$str .= "\t";
		$str .= "$gene";
		$str .= "\t$type\t$aaref\t$aaalt\t$prediction\t";
		#warn $hpocounts[0],"\n";
		$str .= "$hpocounts[0]\t$hpocounts[1]\t$hpocounts[2]\t";
		@pm = &uniq(@pm);
		foreach(@pm){
			$str .= "$_,";
		}
		chop($str) if $str =~ /,$/;
		$str .= "\n";
		push (@report, $str);
		#warn "one loop\n";
	}
	close(VF);
}

#call ANNOVAR scripts
sub annovar()
{
	my $vcfin = shift;
	my $out = shift;
	my $annovar_home = "$Bin/../annovar";
	if ($annovar_home eq '\n'){
		print "Please Set \$ANNOVAR_HOME to your annovar home path\n";
	}
	$annovar_home =~ s/\n//;
	`perl $annovar_home/convert2annovar.pl -format vcf4 $vcfin > $avi_file`;
	`perl $annovar_home/annotate_variation.pl -out $out -build hg19 -hgvs $avi_file $annovar_home/humandb/`;
}

sub filterfamily
{
	my $vcfin = shift;
	my $patientstr = shift;
	my $normalstr = shift;
	my @patients = split(/,/,$patientstr);
	my @normals = split(/,/,$normalstr);
	my @patientcols;
	my @normalcols;
	open (OUT, '>', $ffin) or die "Can't open $ffin:$!\n";
	open (IN, '<', $vcfin) or die "Cant't open $vcfin:$!\n";
	while(<IN>){
		chomp;
		if (/^#CHROM/i){
			my @head = split(/\t/);
			for (my $i =0;$i <=$#head;$i++){
				foreach (@patients){
					if ($head[$i] =~ /$_/i){
						push(@patientcols,$i);
					}
				}
				foreach (@normals){
					if ($head[$i] =~ /$_/i){
						push(@normalcols,$i);
					}
				}
			}
			print OUT join("\t",@head),"\n";
		}elsif(/^#/){
			print OUT "$_\n";
		}else{
			my @line = split(/\t/);
			my $firstcol = $patientcols[0];
			my $firstgp = $line[$firstcol];
			$firstgp =~ s/:.*$//g;
			foreach(@patientcols){
				my $gp = $line[$_];
				$gp =~ s/:.*$//g;
				if ($gp ne $firstgp){
					goto NEXT;
				}
			}
		 	if ($firstgp =~ /^([1-9]\/[1-9])/){
				foreach(@normalcols){
					my $genotype = $1 if ($line[$_] =~ /^([0-9]\/[0-9]):/);
					if (defined $genotype and $genotype eq $firstgp){
						goto NEXT;
					}
				}
				print OUT join("\t",@line),"\n";
			}else{
				next;
			}
		}
		NEXT:
	}
	close(IN);
	close(OUT); 
}

#filter high maf variants
sub filter()
{
	my $vcfin = shift;
	my $maf_rate = shift;
	my @vcf;
	open(MAF,'<',$MAF_FILE);
	my $str = do{local $/=undef;<MAF>;};
	$str =~ s/\n/\t/g;
	my @mlines = split(/\t/,$str);
	my %hash = (@mlines);
	close(MAF);
	open(FILTER, '>', $filter_vcf) or die "Can't open $filter_vcf:$!";
	my $total = `wc -l $vcfin|cut -d ' ' -f 1`;
	chomp($total);
	open(VCF, '<', $vcfin) or die "Can't open $vcfin:$!";
	while(my $line = <VCF>){
		my $flag = 0;
		$percent += 1/$total * 0.1;
		printf LOG "%d%%\n",$percent * 100;
		if ($line =~ /^#/){
			print FILTER $line;
			next;
		}else{
			$htmlhash{"variants_all"} += 1;
			chomp($line);
			my @cols = split(/\t/,$line);
			my $chr = $cols[0];
			$chr =~ s/chr// if $chr =~ /chr/;
			if( $cols[4] =~ /,/){
				my @alts = split(/,/,$cols[4]);
				my @palts;
				for my $alt (@alts){
					my $key;
					if(length($alt) == length($cols[3]) and substr($alt,1,-1) eq substr($cols[3],1,-1)){
						my $ref = substr($cols[3], 0,1);
						my $nalt = substr($alt, 0,1);
						$key = "$chr\_$cols[1]\_$ref\_$nalt";
					}else{
						$key = "$chr\_$cols[1]\_$cols[3]_$alt";
					}
					if (exists $hash{$key}){
						if($hash{$key} > $maf_rate){
							next;
						}else{
							push (@palts, $alt);
						}
					}else{
						push(@palts, $alt);
					}
				}
				if(@palts){
					
					$cols[4] = join(",",@palts);
					$line = join("\t",@cols);
					print FILTER "$line\n";
				}
			}else{
				my $alt = $cols[4];
				my $key = "$chr\_$cols[1]\_$cols[3]_$alt";
				if (exists $hash{$key}){
					if($hash{$key} > $maf_rate){
						next;
					}else{
						print FILTER "$line\n";
				
					}
				}else{
					print FILTER "$line\n";
				}
			}
		}
	}
	close(VCF);
	close(FILTER);
}


sub exportcsv()
{
	my $out = shift;
	my $rank = 20;
	my $outhgmd =  $out . "/html/data/hgmd.csv";
	my $outmed = $out . "/html/data/med.csv";
	my $outhpo = $out . "/html/data/hpo.csv";
	my $outwords = $out . "/html/data/bubble.csv";
	open(OUT1, '>', $outhgmd) or die "Can't open $outhgmd:$!";
	open(OUT2, '>', $outmed) or die "Can't open $outmed:$!";
	open(OUT3, '>', $outhpo) or die "Can't open $outhpo:$!";
	print OUT1 "HPOTERM,DISEASE,GENE,TYPE,PREDICTION,HPO\n";
	print OUT2 "HPOTERM,DISEASE,GENE,TYPE,PREDICTION,HPO\n";
	print OUT3 "HPOTERM,DISEASE,GENE,TYPE,PREDICTION,HPO\n";
	for my $key (keys %hpohash){
		my $value = $hpohash{$key};
		my @ele = split(/\t/, $key);
		next if $ele[-1] =~ /Unknown/;
		next if $ele[-1] =~ /benign/;
		$key =~ s/,/;/g;
		$key =~ s/\t/,/g;
		if ($$value{HGMD} != 0){
			print OUT1 $key;
			print OUT1 ",$$value{HGMD}\n";
		}elsif($$value{GENE_MED} != 0){
			print OUT2 $key;
			print OUT2 ",$$value{GENE_MED}\n";
		}elsif($$value{GENE_HPO} != 0){
			print OUT3 $key;
			print OUT3 ",$$value{GENE_HPO}\n";
		}
	}
	close(OUT1);
	close(OUT2);
	close(OUT3);
	open(OUT, '>', $outwords) or die "Can;t open $outwords:$!";
	print OUT "keyword,count\n";
	for my $key (sort {$keywords{$b} <=> $keywords{$a}}  keys %keywords){
		last if ($rank == 0);
		print OUT "$key,$keywords{$key}\n";
		$rank -= 1;
	}
	close(OUT);
}
sub exporttable()
{
	my $out = shift;
	my $outtable = $out . "/table.txt";
	my @hgmd;
	my @med;
	my @hpo;
	my @other;
	my @lines;
	for my $line (@report){
		my @ele = split(/\t/, $line);
		if ($ele[10] != 0){
			push (@hgmd, $line);
		}elsif ($ele[11] != 0){
			push (@med, $line);
		}elsif ($ele[12] != 0){
			push (@hpo, $line);
		}else{
			push (@other, $line);
		}
	}
	push (@lines, @hgmd);
	push (@lines, @med);
	push (@lines, @hpo);
	push (@lines, @other);
	open (OUT, '>', $outtable) or die "Can't open $outtable:$!";
	print OUT "#chr\tpos\tref\talt\tdisease\tgene\ttype\taa1\taa2\tprediction\tHGMD\tMED\tHPO\tPUBMED\n";
	for my $line (@lines){
		$htmlhash{"variants_anno"} += 1;
		my @ele = split(/\t/,$line);
		my $chr = $ele[0];
		if($chr eq "chr1" or $chr eq "1"){
			$htmlhash{"var_chr1_an"} += 1;
		}elsif($chr eq "chr2" or $chr eq "2"){
			$htmlhash{"var_chr2_an"} += 1;
		}elsif($chr eq "chr3" or $chr eq "3"){
			$htmlhash{"var_chr3_an"} += 1;
		}elsif($chr eq "chr4" or $chr eq "4"){
			$htmlhash{"var_chr4_an"} += 1;
		}elsif($chr eq "chr5" or $chr eq "5"){
			$htmlhash{"var_chr5_an"} += 1;
		}elsif($chr eq "chr6" or $chr eq "6"){
			$htmlhash{"var_chr6_an"} += 1;
		}elsif($chr eq "chr7" or $chr eq "7"){
			$htmlhash{"var_chr7_an"} += 1;
		}elsif($chr eq "chr8" or $chr eq "8"){
			$htmlhash{"var_chr8_an"} += 1;
		}elsif($chr eq "chr9" or $chr eq "9"){
			$htmlhash{"var_chr9_an"} += 1;
		}elsif($chr eq "chr10" or $chr eq "10"){
			$htmlhash{"var_chr10_an"} += 1;
		}elsif($chr eq "chr11" or $chr eq "11"){
			$htmlhash{"var_chr11_an"} += 1;
		}elsif($chr eq "chr12" or $chr eq "12"){
			$htmlhash{"var_chr12_an"} += 1;
		}elsif($chr eq "chr13" or $chr eq "13"){
			$htmlhash{"var_chr13_an"} += 1;
		}elsif($chr eq "chr14" or $chr eq "14"){
			$htmlhash{"var_chr14_an"} += 1;
		}elsif($chr eq "chr15" or $chr eq "15"){
			$htmlhash{"var_chr15_an"} += 1;
		}elsif($chr eq "chr16" or $chr eq "16"){
			$htmlhash{"var_chr16_an"} += 1;
		}elsif($chr eq "chr17" or $chr eq "17"){
			$htmlhash{"var_chr17_an"} += 1;
		}elsif($chr eq "chr18" or $chr eq "18"){
			$htmlhash{"var_chr18_an"} += 1;
		}elsif($chr eq "chr19" or $chr eq "19"){
			$htmlhash{"var_chr19_an"} += 1;
		}elsif($chr eq "chr20" or $chr eq "21"){
			$htmlhash{"var_chr20_an"} += 1;
		}elsif($chr eq "chr21" or $chr eq "21"){
			$htmlhash{"var_chr21_an"} += 1;
		}elsif($chr eq "chr22" or $chr eq "22"){
			$htmlhash{"var_chr22_an"} += 1;
		}elsif($chr eq "chrX" or $chr eq "X"){
			$htmlhash{"var_chrX_an"} += 1;
		}elsif($chr eq "chrY" or $chr eq "Y"){
			$htmlhash{"var_chrY_an"} += 1;
		}
		print OUT $line;
	}
	close(OUT);
}

sub exportjson
{
	my $out = shift;
	my %jhash;
	my %hash;
	for my $key(keys %hpohash){
		my $value = $hpohash{$key};
		my @ele = split(/\t/, $key);
		next if $ele[-1] =~ /Unknown/;
		next if $ele[-1] =~ /benign/;
		my $hpo = $ele[0];
		my $disease = $ele[1];
		my $gene = $ele[2];
		my $count = $hpohash{$key}{"HGMD"} + $hpohash{$key}{"GENE_MED"} + $hpohash{$key}{"GENE_HPO"};
		if (exists $hash{$hpo}){
			if (exists $hash{$hpo}{$disease}){
				if (exists $hash{$hpo}{$disease}{$gene}){
					$hash{$hpo}{$disease}{$gene} += $count;
				}else{
					$hash{$hpo}{$disease}{$gene} = $count;
				}
			}else{
				$hash{$hpo}{$disease}{$gene} = $count;
			}
		}else{
			$hash{$hpo}{$disease}{$gene} = $count;
		}
	}
	$jhash{name} = "hpo";
	my @childrentop;
	for my $key(keys %hash){
		my @children;
		for my $key2 (keys %{$hash{$key}}){
			my @children2;
			for my $key3 (keys %{$hash{$key}{$key2}}){
				push (@children2, {name=>$key3, size=>$hash{$key}{$key2}{$key3}});
			}
			push (@children, {name => $key2, children => [@children2]});
		}
		push (@childrentop,{name => $key, children => [@children]});
	}
	$jhash{children} = [@childrentop]; 
	my $json = new JSON;
	my $json_text = $json->encode(\%jhash);
	my $json_file = $out . "/html/data/bar.json";
	open (OUT ,'>', $json_file) or die "Can't open $json_file:$!\n";
	print OUT $json_text;
	close(OUT);
	my @table;
	for my $line(@report){
		chomp($line);
		my @ele = split(/\t/,$line);
		push (@table ,{
			'chr' => $ele[0],
			'pos' => $ele[1],
			'ref' => $ele[2],
			'alt' => $ele[3],
			'disease' => $ele[4],
			'gene' => $ele[5],
			'type' => $ele[6],
			'aa1' => $ele[7],
			'aa2' => $ele[8],
			'prediction' => $ele[9],
			'HGMD' => $ele[10],
			'MED' => $ele[11],
			'HPO' => $ele[12],
			'PUBMED' => $ele[13]
		});
	}
	$json_text = $json->encode(\@table);
	my $table_json = $out . "/html/data/table.json";
	open (OUT ,'>', $table_json) or die "Can't open $table_json:$!\n";
	print OUT $json_text;
	close(OUT);
}

sub exportpara
{
	my $out = shift;
	my $outpara =  $out . "/html/data/parallel.csv";
	open(OUT,'>',$outpara) or die "Can't open $outpara:$!";
	print OUT "Source,Gene,Disease,Phenotype\n";
	for my $key(keys %hpohash){
		my $value = $hpohash{$key};
		$key =~ s/,/;/g;
		my @ele = split(/\t/, $key);
		my $hpoterm = $ele[0];
		my $disease = $ele[1];
		my $gene = $ele[2];
		my $source;
		my $count;
		if ($$value{HGMD} != 0){
			$source = "HGMD";
			$count = $$value{HGMD};
		}elsif($$value{GENE_MED} != 0){
			$source = "MedGen";
			$count = $$value{GENE_MED};
		}elsif($$value{GENE_HPO} != 0){
			$source = "HPO";
			$count = $$value{GENE_HPO};
		}
		while($count > 0){
			print OUT "$source,$gene,$disease,$hpoterm\n";
			$count--;
		}
	}
	close(OUT);
}

sub check
{
	my $vcfin = shift;
	my $err = '';
	my $warn = '';
	open(IN,'<',$vcfin) or die "Can't open $vcfin:$!";
	while(<IN>){
		if(/^#/){
			if (/refernce/){
				unless (/hg19/){
					$err .= "Not validated reference\n";
				}
			}
		}else{
			chomp;
			my @line = split(/\t/,$_);
			if ($line[6] =~ /\./){
				$warn .= "No VQSR\n";
			}
			last;
		}
	}
	close(IN);
	return ($err, $warn);
}

#main function
sub GPSanno()
{
	printf LOG  "%d%%\n","0";
	my ($vcfin, $out, $maf, $patient, $normal) = &getopts(@ARGV);
	my ($err,$warn) = &check($vcfin);
	if ($err){
		open(ERR,'>',"error.log");
		print ERR $err;
		close(ERR);
		exit 1;
	}
	if($warn){
		open(WARN, '>', "warn.log");
		print WARN $warn;
		close(WARN);
	}
	my $cmdline = "./GPSanno -i $vcfin -o $out -maf $maf";
	$htmlhash{"cmd"} = $cmdline;
	my ($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst) = localtime;
	$year += 1990;
	$mon += 1;
	my $date = "$year-$mon-$day $hour:$min";
	$htmlhash{"data"} = $date;
	if ($patient and $normal){
		&filterfamily($vcfin, $patient, $normal);
		&filter($ffin ,$maf);
	}else{
		&filter($vcfin, $maf);
	}
	$percent = sprintf("%.1f",$percent);
	print LOG ($percent * 100),"\n";
	&annovar($filter_vcf, $annovarout);
	$casehash = &loadCASE($table);
	&anno($annovarout);
	$percent = sprintf("%.1f",$percent);
	print LOG ($percent * 100) ,"\n";
	`cp -r $Bin/web $out`;
	&exporttable($out);
	$percent += 1/50;
	printf LOG "%d%%\n",$percent * 100;
	&exportcsv($out);
	$percent += 1/50;
	printf LOG "%d%%\n",$percent * 100;
	&exporthtml($out, \%htmlhash);
	$percent += 1/50;
	printf LOG "%d%%\n",$percent * 100;
	&exportjson($out);
	$percent += 1/50;
	printf LOG "%d%%\n",$percent * 100;
	&exportpara($out);
	$percent += 1/50;
	printf LOG "%d%%\n",$percent * 100;
	`rm temp.*`;
	print "\n";
}
&GPSanno();
close (LOG);
