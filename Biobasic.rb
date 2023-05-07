#Practice#
#Bioinformatics Basic#

File.open("file.txt") do |f|
    while line = f.gets
    print line
    end
    end

File.open("seq.fa") do |in|
    File.open("seq2.fa", "w+") do |out|
    while line = in.gets
    out << line unless line =~
    /^>/
    end
    end
    end

gene_counts = ("Human" => 31000, "Fruit fly" => 13000,
        "Mouse" => 30000, "Chickenpox virus" => 69, "Rice" => 40000,
        "Tuberculosis bacteria" => 4000)
gene_counts.each_pair {|key, value| puts "#{key}has #{value}
        genes in its genome."}

@dict = {"A" => "Adenine", "T" => "Thymine", "G" => "Guanine",
        "C" => "Cytosine"}
sequence = "CTATGCGGTA"
sequence.scan(/./).each {|i| puts @dict[i]}

@sequence="ATGAATCCAAGCCAAATACTTGAAAATTTAAAAAAAGAATTAAGTGAAAAC
GAATACGAAAACTATTTATCAAATTTAAAATTCAACGAAAAACAAAGCAAAGCAGATCTTTT
AGTTTTTAATGCTCCAAATGAACTCATGGCTAAATTCATACAAACAAAATACGGCAAAAAAA
TCGCGCATTTTTATGAAGTGCAAAGCGGAAATAAAGCCATCATAAATATACAAGCACAAAGT
GCTAAACAAAGCAACAAAAGCACAAAAATCGACATAGCTCATATAAAAGCACAAAGCACG";
@sum=0
@sequence.each_char {|i| @sum+=1 if i == 'A'}
puts @sum

@dict = {"TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" =>
"L", "CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
"ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M", "GTT"
=> "V", "GTC" => "V", "GTA" => "V", "GTG" => "V", "TCT" =>
"S", "TCC" => "S", "TCA" => "S", "TCG" => "S", "CCT" => "P",
"CCC" => "P", "CCA" => "P", "CCG" => "P", "ACT" => "T", "ACC"
=> "T", "ACA" => "T", "ACG" => "T", "GCT" => "A", "GCC" =>
"A", "GCA" => "A", "GCG" => "A", "TAT" => "Y", "TAC" => "Y",
152
Ruby vs. Perl – the Languages of Bioinformatics
"TAA" => "STOP", "TAG" => "STOP", "CAT" => "H", "CAC" => "H",
"CAA" => "Q", "CAG" => "Q", "AAT" => "N", "AAC" => "N", "AAA"
=> "K", "AAG" => "K", "GAT" => "D", "GAC" => "D", "GAA" =>
"E", "GAG" => "E", "TGT" => "C", "TGC" => "C", "TGA" =>
"STOP", "TGG" => "W", "CGT" => "R", "CGC" => "R", "CGA" =>
"R", "CGG" => "R", "AGT" => "S", "AGC" => "S", "AGA" => "R",
"AGG" => "R", "GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG"
=> "G"};
@sequence="atgagttctgactctgagatggccatttttggggaggctgctcctttcct
ccgaaagtctgaaagggagcgaattgaagcccagaacaagccttttgatgccaagacatca
gtctttgtggtggaccctaaggagtcctttgtgaaagcaacagtgcagagcagggaagggg
ggaaggtgacagctaagaccgaagctggagctactgtaacagtgaaagatgaccaagtctt
ccccatgaaccctcccaaatatgacaagatcgaggacatggccatgatgactcatctacac
gagcctgctgtgctgtacaacctcaaagagcgctacgcagcctggatgatctacacctact
caggc";
@goal="MSSDSEMAIFGEAAPFLRKSERERIEAQNKPFDAKTSVFVVDPKESFVKATVQS
REGGKVTAKTEAGATVTVKDDQVFPMNPPKYDKIEDMAMMTHLHEPAVLYNLKERYAAWMI
YTYSG";
@translation=""
@sequence.upcase.scan(/.../).each
{
|i| @translation << @dict[i]
}
puts "Success!" if @translation == @goal

