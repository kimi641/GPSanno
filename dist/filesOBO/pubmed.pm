
#!/usr/bin/perl -w
use LWP::Simple;
$query = 'xanthomatosis';
$base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
$url = $base . "esearch.fcgi?db=pubmed&term=$query&usehistory=y";
#post the esearch URL
$output = get($url);
#parse WebEnv, QueryKey and Count (# records retrieved)
$web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
$key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
$count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);
#open output file for writing
open (SALIDA, ">pubmedCTX.txt") || die "Can't open file!\n";
#retrieve data in batches of 500
$retmax = 1000;
for ($retstart = 0; $retstart < $count; $retstart += $retmax) {
        $efetch_url = $base ."efetch.fcgi?db=pubmed&WebEnv=$web";
        $efetch_url .= "&query_key=$key&retstart=$retstart";
        $efetch_url .= "&retmax=$retmax&rettype=medline&retmode=text";
        $efetch_out = get($efetch_url);
        print SALIDA "$efetch_out";
}
close SALIDA;