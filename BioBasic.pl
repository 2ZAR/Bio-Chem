#Practice#
#Bioinformatics Basic#

open(F, "file.txt");
while($line = <F>)
{
print "$line"
}
close F;

open (F, "seq.fa");
open (FF,">seq2.txt");
while (<F>)
{
next if(/^>/);
print FF;
}
close F;
close FF;

%gene_counts = ("Human" => 31000, "Fruit fly" => 13000,
"Mouse" => 30000, "Chickenpox virus" => 69, "Rice" => 40000,
"Tuberculosis bacteria" => 4000);
while ( ( $key, $value ) = each %gene_counts )
{
print "$key has $value genes in its genome.\n";
}

%dict = (A => Adenine, T => Thymine, G => Guanine,
C => Cytosine);
$sequence = 'CTATGCGGTA';
while ( $sequence =~
/./g )
{
print "$dict{$&}\n";
}

$sequence="ATGAATCCAAGCCAAATACTTGAAAATTTAAAAAAAGAATTAAGTGAAAAC
GAATACGAAAACTATTTATCAAATTTAAAATTCAACGAAAAACAAAGCAAAGCAGATCTTTT
AGTTTTTAATGCTCCAAATGAACTCATGGCTAAATTCATACAAACAAAATACGGCAAAAAAA
TCGCGCATTTTTATGAAGTGCAAAGCGGAAATAAAGCCATCATAAATATACAAGCACAAAGT
GCTAAACAAAGCAACAAAAGCACAAAAATCGACATAGCTCATATAAAAGCACAAAGCACG";
$sum=0;
@tab=split('', $sequence);
foreach $i (@tab)
{
$sum++ if $i eq 'A';
}
print $sum;

%dict = ("TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" =>
"L", "CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
151
Maciej Goliński and Agnieszka Kitlas Golińska
"ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M", "GTT"
=> "V", "GTC" => "V", "GTA" => "V", "GTG" => "V", "TCT" =>
"S", "TCC" => "S", "TCA" => "S", "TCG" => "S", "CCT" => "P",
"CCC" => "P", "CCA" => "P", "CCG" => "P", "ACT" => "T", "ACC"
=> "T", "ACA" => "T", "ACG" => "T", "GCT" => "A", "GCC" =>
"A", "GCA" => "A", "GCG" => "A", "TAT" => "Y", "TAC" => "Y",
"TAA" => "STOP", "TAG" => "STOP", "CAT" => "H", "CAC" => "H",
"CAA" => "Q", "CAG" => "Q", "AAT" => "N", "AAC" => "N", "AAA"
=> "K", "AAG" => "K", "GAT" => "D", "GAC" => "D", "GAA" =>
"E", "GAG" => "E", "TGT" => "C", "TGC" => "C", "TGA" =>
"STOP", "TGG" => "W", "CGT" => "R", "CGC" => "R", "CGA" =>
"R", "CGG" => "R", "AGT" => "S", "AGC" => "S", "AGA" => "R",
"AGG" => "R", "GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG"
=> "G");
$sequence=uc("atgagttctgactctgagatggccatttttggggaggctgctccttt
cctccgaaagtctgaaagggagcgaattgaagcccagaacaagccttttgatgccaagaca
tcagtctttgtggtggaccctaaggagtcctttgtgaaagcaacagtgcagagcagggaag
gggggaaggtgacagctaagaccgaagctggagctactgtaacagtgaaagatgaccaagt
cttccccatgaaccctcccaaatatgacaagatcgaggacatggccatgatgactcatcta
cacgagcctgctgtgctgtacaacctcaaagagcgctacgcagcctggatgatctacacct
actcaggc");
$goal="MSSDSEMAIFGEAAPFLRKSERERIEAQNKPFDAKTSVFVVDPKESFVKATVQS
REGGKVTAKTEAGATVTVKDDQVFPMNPPKYDKIEDMAMMTHLHEPAVLYNLKERYAAWMI
YTYSG";
$translation="";
while ($sequence =~
/.../g)
{
$translation .= $dict{$&};
}
print "Success!" if ($translation eq $goal);


