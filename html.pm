#!/usr/bin/perl -w
use strict;
sub exporthtml
{
	my $out = shift;
	my $phash = shift;
	chop($out) if ($out =~ /\/$/);
	my $html = $out . "/summary.html";
	open (OUT, '>', $html) or die "Can't open $html:$!";
	print OUT "<head>\n";
	print OUT "<link rel=\"stylesheet\" href=\"./lib/styles/summary.css\">\n";
	print OUT "</head>\n";
	print OUT "<center>\n";	
	print OUT "<div class=\"main\">\n";
	print OUT "<center><h3>GPSanno Summary</h3></center>\n";
	print OUT "<center><div style=\"margin-left: .5em\">\n";
	print OUT "<table class=\"toc\"><tr><td>\n";
	print OUT "<center><b>Contents</b></center>\n";
	print OUT "<a href=\"#summary\">Summary</a><br>\n";
	print OUT "<a href=\"#varChr\">Variants detail by chromsome</a><br>\n";
	print OUT "<a href=\"#varType\">Variants by type</a><br>\n";
	print OUT "<a href=\"#varPrediction\">Number of variants by prediciton</a><br>\n";
	print OUT "<a href=\"html/Weightedtree.html\">Weightedtree for HPO</a><br>\n";
	print OUT "<a href=\"html/bar-hierarchy.html\">Bar-hierarchy for HPO</a><br>\n";
	print OUT "<a href=\"html/parallel.html\">Parallel Set for HPO</a><br>\n";
	#print OUT "<a href=\"html/Bubble.html\">Key word rank</a><br>\n";
	print OUT "<a href=\"html/table.html\">Table report</a><br>\n";
	print OUT "</tr></td></table>\n";
	print OUT "</div></center>\n";
	print OUT "<hr><a name=\"summary\">\n";
	print OUT "<center><b> Summary </b><p>\n";
	print OUT "<table class=\"summary\" border=0>\n";
	print OUT "<tr bgcolor=ffffff><td valign=top><b> Genome </b></td><td> Hg19 </td></tr>\n";
	print OUT "<tr bgcolor=dddddd><td vailgn=top><b> Data </b></td><td> $$phash{data} </td></tr>\n";
	print OUT "<tr bgcolor=ffffff><td valign=top><b> Command line arguments </b></td><td> $$phash{cmd} </td></tr>\n";
	print OUT "<tr bgcolor=dddddd><td vailgn=top><b> Number of variants (all) </b></td><td> $$phash{variants_all} </td></tr>\n";
	print OUT "<tr bgcolor=ffffff><td valign=top><b> Number of variants (keep) </b></td><td> $$phash{variants_fl} </td></tr>\n";
	print OUT "<tr bgcolor=dddddd><td vailgn=top><b> Number of variants (annotated) </b></td><td> $$phash{variants_anno} </td></tr>\n";
	print OUT "</table><p></center>\n";
	print OUT "<hr><a name=\"varChr\">\n";
	print OUT "<center><b> Variants detail by chromsome </b><p>\n";
	print OUT "<table border=1>\n";
	print OUT "<thead><th> Chromosome </th><th> Annotated variants </th><th> Variants </th></thead>\n";
	print OUT "<tbody>\n";
	print OUT "<tr><td> 1 </td><td class=\"numeric\"> $$phash{var_chr1_an} </td><td class=\"numeric\"> $$phash{var_chr1} </td></tr>\n";
	print OUT "<tr><td> 2 </td><td class=\"numeric\"> $$phash{var_chr2_an} </td><td class=\"numeric\"> $$phash{var_chr2} </td></tr>\n";
	print OUT "<tr><td> 3 </td><td class=\"numeric\"> $$phash{var_chr3_an} </td><td class=\"numeric\"> $$phash{var_chr3} </td></tr>\n";
	print OUT "<tr><td> 4 </td><td class=\"numeric\"> $$phash{var_chr4_an} </td><td class=\"numeric\"> $$phash{var_chr4} </td></tr>\n";
	print OUT "<tr><td> 5 </td><td class=\"numeric\"> $$phash{var_chr5_an} </td><td class=\"numeric\"> $$phash{var_chr5} </td></tr>\n";
	print OUT "<tr><td> 6 </td><td class=\"numeric\"> $$phash{var_chr6_an} </td><td class=\"numeric\"> $$phash{var_chr6} </td></tr>\n";
	print OUT "<tr><td> 7 </td><td class=\"numeric\"> $$phash{var_chr7_an} </td><td class=\"numeric\"> $$phash{var_chr7} </td></tr>\n";
	print OUT "<tr><td> 8 </td><td class=\"numeric\"> $$phash{var_chr8_an} </td><td class=\"numeric\"> $$phash{var_chr8} </td></tr>\n";
	print OUT "<tr><td> 9 </td><td class=\"numeric\"> $$phash{var_chr9_an} </td><td class=\"numeric\"> $$phash{var_chr9} </td></tr>\n";
	print OUT "<tr><td> 10 </td><td class=\"numeric\"> $$phash{var_chr10_an} </td><td class=\"numeric\"> $$phash{var_chr10} </td></tr>\n";
	print OUT "<tr><td> 11 </td><td class=\"numeric\"> $$phash{var_chr11_an} </td><td class=\"numeric\"> $$phash{var_chr11} </td></tr>\n";
	print OUT "<tr><td> 12 </td><td class=\"numeric\"> $$phash{var_chr12_an} </td><td class=\"numeric\"> $$phash{var_chr12} </td></tr>\n";
	print OUT "<tr><td> 13 </td><td class=\"numeric\"> $$phash{var_chr13_an} </td><td class=\"numeric\"> $$phash{var_chr13} </td></tr>\n";
	print OUT "<tr><td> 14 </td><td class=\"numeric\"> $$phash{var_chr14_an} </td><td class=\"numeric\"> $$phash{var_chr14} </td></tr>\n";
	print OUT "<tr><td> 15 </td><td class=\"numeric\"> $$phash{var_chr15_an} </td><td class=\"numeric\"> $$phash{var_chr15} </td></tr>\n";
	print OUT "<tr><td> 16 </td><td class=\"numeric\"> $$phash{var_chr16_an} </td><td class=\"numeric\"> $$phash{var_chr16} </td></tr>\n";
	print OUT "<tr><td> 17 </td><td class=\"numeric\"> $$phash{var_chr17_an} </td><td class=\"numeric\"> $$phash{var_chr17} </td></tr>\n";
	print OUT "<tr><td> 18 </td><td class=\"numeric\"> $$phash{var_chr18_an} </td><td class=\"numeric\"> $$phash{var_chr18} </td></tr>\n";
	print OUT "<tr><td> 19 </td><td class=\"numeric\"> $$phash{var_chr19_an} </td><td class=\"numeric\"> $$phash{var_chr19} </td></tr>\n";
	print OUT "<tr><td> 20 </td><td class=\"numeric\"> $$phash{var_chr20_an} </td><td class=\"numeric\"> $$phash{var_chr20} </td></tr>\n";
	print OUT "<tr><td> 21 </td><td class=\"numeric\"> $$phash{var_chr21_an} </td><td class=\"numeric\"> $$phash{var_chr21} </td></tr>\n";
	print OUT "<tr><td> 22 </td><td class=\"numeric\"> $$phash{var_chr22_an} </td><td class=\"numeric\"> $$phash{var_chr22} </td></tr>\n";
	print OUT "<tr><td> X </td><td class=\"numeric\"> $$phash{var_chrX_an} </td><td class=\"numeric\"> $$phash{var_chrX} </td></tr>\n";
	print OUT "<tr><td> Y </td><td class=\"numeric\"> $$phash{var_chrY_an} </td><td class=\"numeric\"> $$phash{var_chrY} </td></tr>\n";
	print OUT "</tbody>\n";
	print OUT "</table><center>\n";
	print OUT "<hr><a name=\"varType\"><center><b> Number of variants by type </b><p>\n";
	print OUT "<table border=1>\n";
	print OUT "<thead><th><b> Type </b></th><th><b> Total </b></th></tr></thead>\n";
	print OUT "<tbody>\n";
	print OUT "<tr><td><b></b> frameshift deletion </td><td class=\"numeric\" bgcolor=\"ffffff\"> $$phash{frameshift_deletion} </td></tr>\n";
	print OUT "<tr><td><b></b> frameshift insertion </td><td class=\"numeric\" bgcolor=\"ffffff\"> $$phash{frameshift_insertion} </td></tr>\n";
	print OUT "<tr><td><b></b> nonframeshift deletion </td><td class=\"numeric\" bgcolor=\"ffffff\"> $$phash{nonframeshift_deletion} </td></tr>\n";
	print OUT "<tr><td><b></b> nonframeshift insertion </td><td class=\"numeric\" bgcolor=\"ffffff\"> $$phash{nonframeshift_insertion} </td></tr>\n";
	print OUT "<tr><td><b></b> nonsynonymous SNV </td><td class=\"numeric\" bgcolor=\"ffffff\"> $$phash{nonsynonymous_SNV} </td></tr>\n";
	print OUT "<tr><td><b></b> stopgain </td><td class=\"numeric\" bgcolor=\"ffffff\"> $$phash{stopgain} </td></tr>\n";
	print OUT "<tr><td><b></b> stoploss </td><td class=\"numeric\" bgcolor=\"ffffff\"> $$phash{stoploss} </td></tr>\n";
	print OUT "<tr><td><b></b> synonymous SNV </td><td class=\"numeric\" bgcolor=\"ffffff\"> $$phash{synonymous_SNV} </td></tr>\n";
	print OUT "<tr><td><b></b> unknown </td><td class=\"numeric\" bgcolor=\"ffffff\"> $$phash{unknown} </td></tr>\n";
	print OUT "</tbody>\n";
	print OUT "</table></center>\n";
	print OUT "<hr><a name=\"varPrediction\"><center><b> Number of variants by prediction </b><p>\n";
	print OUT "<table border=1>\n";
	print OUT "<thead><tr><th><b> Prediction </b></th><th> Count </th></tr></thead>\n";
	print OUT "<tr><td><b> probably damaging </b></td><td class=\"numeric\" bgcolor=\"ffffff\"> $$phash{probably_damaging} </td></tr>\n";
	print OUT "<tr><td><b> possibly damaging </b></td><td class=\"numeric\" bgcolor=\"ffffff\"> $$phash{possibly_damaging} </td></tr>\n";
	print OUT "<tr><td><b> benign </b></td><td class=\"numeric\" bgcolor=\"ffffff\"> $$phash{benign} </td></tr>\n";
	print OUT "</table><br><p></center>\n";
	print OUT "</div></center>\n";
	close(OUT);
}
1;

