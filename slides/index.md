
class: center, middle






# Processing sequence format file using BioJulia






### Taipei Bioinformatics Omnibus






#### Yueh-Hua Tu






#### 2020.10.24


---






# Outline


  * Processing fastq/fasta format files for DNA/protein
  * Processing sam/bam format files for alignment
  * pipeline and shell


---






# Constructing DNA/RNA/amino acid sequence


```julia
using BioSequences
```






### Construction from strings


  * `dna"ATCG"`
  * `LongDNASeq("TTANC")`
  * `LongSequence{DNAAlphabet{4}}("TTANC")`


--


  * `rna"AUGC"`
  * `LongRNASeq("UUANC")`


--


  * `aa"ARNDCQEGHI"`
  * `LongAminoAcidSeq("ARNDCQEGHI")`


--


  * `ReferenceSequence("NNCGTATTTTCN")`


---




# Constructing DNA/RNA/amino acid sequence






### Construction from BioSymbols


  * `LongDNASeq([DNA_T, DNA_T, DNA_A, DNA_N, DNA_C])`
  * `LongRNASeq([RNA_U, RNA_U, RNA_A, RNA_N, RNA_C])`
  * `LongAminoAcidSeq([AA_K, AA_I, AA_A, AA_C, AA_L])`


--






### Concatenation with sequences


  * `LongDNASeq(dna"ATCG", dna"NNN", dna"TCGA")`
  * `dna"ATCG" * dna"TCGA"`
  * `repeat(dna"ATCG", 10)`
  * `dna"ATCG"^10`


---






# Long Sequences


`LongSequence` is an efficient, general purpose sequence type for creating and editing. `LongSequence` is composed of `Alphabet` and acts as a container. Precisely, `LongSequence` is parametrized by a `Alphabet` type.


  * `LongSequence{DNAAlphabet{4}}`

      * alias: `LongDNASeq`
  * `LongSequence{RNAAlphabet{4}}`

      * alias: `LongRNASeq`
  * `LongSequence{AminoAcidAlphabet}`

      * alias: `LongAminoAcidSeq`
  * `LongSequence{CharAlphabet}`

      * alias: `LongCharSeq`


---






# BioSymbols






### Alphabet


`DNAAlphabet{2}` is an alphabet that stores a base with two bits and it allows only unambiguous nucleotide symbols. `DNAAlphabet{4}` allows sequences with ambiguous nucleotides. It allows you to choose the memory layout as you need.


  * `ACGT` in DNA
  * `ACGU` in RNA


--






### APIs


  * `alphabet(DNA)`
  * `alphabet(RNA)`
  * `alphabet(AminoAcid)`
  * `gap(DNA)`, `gap(RNA)` and `gap(AminoAcid)`


---






# APIs for Alphabet






### `isambiguous`


```julia
julia> isambiguous(DNA_N)
true

julia> isambiguous(DNA_A)
false
```


--






### `iscompatible`


```julia
julia> iscompatible(DNA_A, DNA_A)
true

julia> iscompatible(DNA_C, DNA_N)  # DNA_N can be DNA_C
true

julia> iscompatible(DNA_N, DNA_C)  # Commutative
true

julia> iscompatible(DNA_C, DNA_R)  # DNA_R can be A or G, not C
false
```


---






# Indexing sequences


```julia
julia> seq = dna"ATCG"
4nt DNA Sequence:
ATCG

julia> seq[2]
DNA_T

julia> seq[2:end]
3nt DNA Sequence:
TCG

julia> seq[begin:4]
4nt DNA Sequence:
ATCG
```


--


```julia
julia> length(seq)
4
```


---






# Modifying sequences


```julia
julia> seq[4] = DNA_A
DNA_A

julia> seq
4nt DNA Sequence:
ATCA
```


--






### Modifying subsequence doesn't affect the orginal


```julia
julia> seq = dna"ATATATATATATATATAT"
18nt DNA Sequence:
ATATATATATATATATAT

julia> subseq = seq[1:5]
5nt DNA Sequence:
ATATA
```


--


```julia
julia> subseq[3] = DNA_C
DNA_C

julia> subseq
5nt DNA Sequence:
ATCTA

julia> seq
18nt DNA Sequence:
ATATATATATATATATAT
```


  * *copy-on-write* strategy


---




# Modifying sequences


```julia
julia> push!(seq, DNA_A)
19nt DNA Sequence:
ATATATATATATATATATA

julia> pop!(seq)
DNA_A

julia> seq
18nt DNA Sequence:
ATATATATATATATATAT
```


--


```julia
julia> pushfirst!(seq, DNA_C)
19nt DNA Sequence:
CATATATATATATATATAT

julia> popfirst!(seq)
DNA_C
```


--


```julia
julia> insert!(seq, 2, DNA_G)
19nt DNA Sequence:
AGTATATATATATATATAT
```


---






# Modifying sequences - more APIs


  * `deleteat!`
  * `append!`


--


  * `reverse!`, `reverse`
  * `complement!`, `complement`
  * `reverse_complement!`, `reverse_complement`


--


  * `ungap!`, `ungap`
  * `canonical!`, `canonical`


---






# Transforming sequences






### Transcription


```julia
julia> dna = dna"ATATATATATATATATAT"
18nt DNA Sequence:
ATATATATATATATATAT

julia> rna = convert(LongRNASeq, dna)
18nt RNA Sequence:
AUAUAUAUAUAUAUAUAU
```


--






### Translation


```julia
julia> translate(rna)
6aa Amino Acid Sequence:
IYIYIY
```


---




# Translation






### Genetic code list


```julia
julia> ncbi_trans_table
Translation Tables:
  1. The Standard Code (standard_genetic_code)
  2. The Vertebrate Mitochondrial Code (vertebrate_mitochondrial_genetic_code)
  3. The Yeast Mitochondrial Code (yeast_mitochondrial_genetic_code)
  ...
```


--






### Translation API


Translates `LongRNASeq` or `LongDNASeq` to `LongAminoAcidSeq`.


```julia
translate(seq, code=standard_genetic_code, allow_ambiguous_codons=true, convert_start_codon=false)
```


---






# Checking properties of your sequences


```julia
julia> isrepetitive(dna"AAGT", 2)
true

julia> isrepetitive(dna"AAAGT", 3)
true
```


--


```julia
julia> ispalindromic(dna"GAATTC")
true
```


--


```julia
julia> hasambiguity(dna"ATNGC")
true
```


---






# Generating random sequences


```julia
julia> randseq(DNAAlphabet{2}(), 27)
27nt DNA Sequence:
TGACGTCATTGTTTCCAGACCGGAGAG

julia> randseq(RNAAlphabet{4}(), 27)
27nt RNA Sequence:
ACCUAUCGAAAUCUCUACAUUGCGCUC

julia> randseq(AminoAcidAlphabet(), 50)
50aa Amino Acid Sequence:
FYINRMCKECCWQMKVHSYSYAPNQIDGIIAMATWMISEKQSAYQDHMRL
```


--






### Shortcut


  * `randdnaseq`
  * `randrnaseq`
  * `randaaseq`


---




# Generating random sequences






### Setting sampling weight


Uniformly distributed sampler


```julia
julia> sp = SamplerUniform(dna"ACGT")
SamplerUniform{DNA}(DNA[DNA_A, DNA_C, DNA_G, DNA_T])

julia> seq = randseq(DNAAlphabet{2}(), sp, 100)
100nt DNA Sequence:
ACCTTAACTGACGGGGCAAGAACCTCATCAAGCGGCATA…ACACTCTGACCTGCGTGAGTATCCCTTAAATGCGCGCTA
```


--


Weighted sampler


```julia
julia> sp = SamplerWeighted(dna"ACGTN", [0.1, 0.3, 0.3, 0.3])
SamplerWeighted{DNA}(DNA[DNA_A, DNA_C, DNA_G, DNA_T, DNA_N], [0.1, 0.3, 0.3, 0.3, 0.0])

julia> seq = randseq(DNAAlphabet{2}(), sp, 100)
100nt DNA Sequence:
CCGTAGATTCGGGCTGGCGGGCGTTTACTGCGCTTGGGC…TCCTGCGTTTTGCTGGGGCTGTTTGTGGTGTGTGCTCCC
```


---






# Download sequence data from NCBI


```julia
using BioServices

res = BioServices.EUtils.einfo(db="pubmed")
```


--


```julia
write("pubmed.xml", res.body)
```


--


Or you may want to parse it into a tree structure.


```julia
julia> using EzXML

julia> doc = parsexml(res.body)
EzXML.Document(EzXML.Node(<DOCUMENT_NODE@0x000055b3d1738670>))

julia> root(doc)
EzXML.Node(<ELEMENT_NODE[eInfoResult]@0x000055b3d3234db0>)
```


---






# E-utility in BioServices


  * `einfo`
  * `esearch`
  * `epost`
  * `esummary`
  * `efetch`
  * `elink`


> [Ref: E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25499/)



---






# Download BRCA1 gene sequnece from NCBI


```julia
julia> res = BioServices.EUtils.efetch(db="nuccore", id="MF945608", retmode="xml")
HTTP.Messages.Response:
"""
HTTP/1.1 200 OK
Date: Thu, 22 Oct 2020 02:37:53 GMT
Server: Finatra
...
```


--


```julia
julia> doc = parsexml(res.body)
EzXML.Document(EzXML.Node(<DOCUMENT_NODE@0x000055b3d3062db0>))

julia> seq = findfirst("/GBSet/GBSeq", doc)
EzXML.Node(<ELEMENT_NODE[GBSeq]@0x000055b3d2b20560>)
```


--


```julia
julia> nodecontent(findfirst("GBSeq_definition", seq))
"Homo sapiens isolate 44 BRCA1 (BRCA1) gene, partial cds"
```


--


```julia
julia> nodecontent(findfirst("GBSeq_accession-version", seq))
"MF945608.1"
```


--


```julia
julia> nodecontent(findfirst("GBSeq_sequence", seq))
"ctcactaaagacagaatgaatgtagaaaaggctgaattctgtaataaaagcaaacagcctggcttagcaaggagccaacataacagatgggctggaagtaaggaaacatgtaatgataggcggactcccagcacagaaaaaaaggtagatctgaatgctgatcccctgtgtgagagaaaagaatggaataagcagaaactgccatgctcagagaatcctagagatactgaagatgttccttggataacactaaatagcagcattcagaaag"
```


---






# Get BRCA1 sequence


```julia
julia> using BioSequences

julia> dna = LongDNASeq(nodecontent(findfirst("GBSeq_sequence", seq)))
271nt DNA Sequence:
CTCACTAAAGACAGAATGAATGTAGAAAAGGCTGAATTC…ATGTTCCTTGGATAACACTAAATAGCAGCATTCAGAAAG
```


---






# I/O


```bash
prefetch SRR12876566 -O data/
fastq-dump data/SRR12876566/SRR12876566.sra --split-files SRR12876566
```


> *Salmonella sequenced by MiSeq from SRR12876566.*



--






### FASTQ reader


```julia
using FASTX

filename = "data/SRR12876566/SRR12876566_1.fastq"
r = FASTQ.Reader(open(filename, "r"))
rec = []
for record in r
    push!(rec, record)
end
close(reader)
```


--


```julia
julia> length(rec)
1251960
```


---




### FASTQ reader


```julia
julia> sequence(rec[1])
210nt DNA Sequence:
TCCCTCAAGCGCTCAACTGGACGCCTTATCTTCCCTCAC…TAGTTTCTACTGATAATAATTTTGGGCTTATTATAGCTA
```


--


```julia
julia> quality(rec[1])
210-element Array{UInt8,1}:
 0x10
 0x10
 0x1d
 0x10
 0x1d
```


--






### FASTQ writer


```julia
w = FASTQ.Writer(open("output.fastq", "w"))
write(w, rec)
close(w)
```


--


Or simply...


```julia
w = open(FASTQ.Writer, "output.fastq")
```


---






### FASTA reader


```julia
rec = []
reader = open(FASTA.Reader, "data/salmonella_lt2.fasta")
for record in reader
    push!(rec, record)
end
close(reader)
```


--


Or simply...


```julia
rec = [r for r in open(FASTA.Reader, "data/salmonella_lt2.fasta")]
```


--






#### Convert into DNA sequence


```julia
julia> seq = sequence(LongDNASeq, rec[1])
4857450nt DNA Sequence:
AGAGATTACGTCTGGTTGCAAGAGATCATGACAGGGGGA…AAACTAACAAAATAACGTGCTGTAATTTTTAAAATAATA
```


---






### FASTA writer


```julia
w = FASTA.Writer(open("output.fasta", "w"))
# Or open(FASTA.Writer, "output.fasta")
write(w, records)
close(w)
```


--






### Get information


```julia
julia> identifier(rec[1])
"NC_003197.2"

julia> description(rec[1])
"Salmonella enterica subsp. enterica serovar Typhimurium str. LT2, complete genome"
```


---






# Searching






### Exact matching


```julia
julia> seq
4857450nt DNA Sequence:
AGAGATTACGTCTGGTTGCAAGAGATCATGACAGGGGGA…AAACTAACAAAATAACGTGCTGTAATTTTTAAAATAATA

julia> findfirst(DNA_G, seq)
2

julia> findfirst(dna"TTGCA", seq)
16:20
```


--


```julia
julia> findlast(dna"TTGCA", seq)
4857382:4857386
```


--


```julia
julia> occursin(dna"TTGCA", seq)
true
```


---






### Matching ambiguous mers


```julia
julia> findfirst(dna"CCGA", dna"ACCNA")
2:5
```


--






### Accelerating exact match by preprocessing query sequence


```julia
julia> query = ExactSearchQuery(dna"CCGA")
ExactSearchQuery{LongSequence{DNAAlphabet{4}}}(CCGA, 0x00000007, 4, 1)

julia> findfirst(query, dna"ACCNA")
2:5
```


---






# Approximate search


```julia
julia> approxsearch(seq, dna"TTACG", 1)
6:9
```


--


```julia
julia> approxsearch(seq, dna"TTACG", 2)
6:8
```


--






### Accelerating approximate match by preprocessing query sequence


```julia
julia> query = ApproximateSearchQuery(dna"TTACG")
ApproximateSearchQuery{LongSequence{DNAAlphabet{4}}}(TTACG, UInt64[0x0000000000000000, 0x0000000000000004, 0x0000000000000008, 0x000000000000000c, 0x0000000000000010, 0x0000000000000014, 0x0000000000000018, 0x000000000000001c, 0x0000000000000003, 0x0000000000000007, 0x000000000000000b, 0x000000000000000f, 0x0000000000000013, 0x0000000000000017, 0x000000000000001b, 0x000000000000001f], UInt64[0x0000000000000000, 0x0000000000000004, 0x0000000000000002, 0x0000000000000006, 0x0000000000000001, 0x0000000000000005, 0x0000000000000003, 0x0000000000000007, 0x0000000000000018, 0x000000000000001c, 0x000000000000001a, 0x000000000000001e, 0x0000000000000019, 0x000000000000001d, 0x000000000000001b, 0x000000000000001f], [0, 1, 2, 2, 1, 1])
```


--


```julia
julia> approxsearch(seq, query, 1)
6:9

julia> approxsearch(seq, query, 2)
6:8
```


---






# Pattern matching by regular expression






### BIOlogical REgular expression


```julia
julia> pat = biore"A+"dna
biore"A+"dna
```


--


```julia
julia> match(pat, dna"AAAACC")
RegexMatch("AAAA")
```


--


```julia
julia> m = match(biore"A+C+"dna, dna"AAAACC")
RegexMatch("AAAACC")
```


--






### `occursin`


```julia
julia> occursin(biore"A+C+"dna, dna"AAAACC")
true
```


--


```julia
julia> occursin(biore"A+C+"dna, dna"CCC")
false
```


---






# `RegexMatch` object


```julia
julia> m.seq
6nt DNA Sequence:
AAAACC

julia> m.captured
2-element Array{Int64,1}:
 1
 7
```


--






### `eachmatch`


```julia
julia> collect(matched(x) for x in eachmatch(biore"TATA*?"dna, dna"TATTATAATTA")) # overlap
4-element Array{LongSequence{DNAAlphabet{4}},1}:
 TAT
 TAT
 TATA
 TATAA
```


--


```julia
julia> collect(matched(x) for x in eachmatch(biore"TATA*"d, dna"TATTATAATTA", false)) # no overlap
2-element Array{LongSequence{DNAAlphabet{4}},1}:
 TAT
 TATAA
```


[Ref. Biological regular expression](https://biojulia.net/BioSequences.jl/stable/sequence_search/#Regular-expression-search-1)


---






# Demultiplexing


```julia
julia> barcodes = [dna"ATGG", dna"CAGA", dna"GGAA", dna"TACG"];

julia> dplxr = Demultiplexer(barcodes, n_max_errors=1, distance=:hamming)
Demultiplexer{LongSequence{DNAAlphabet{4}}}:
  distance: hamming
  number of barcodes: 4
  number of correctable errors: 1
```


--


```julia
julia> demultiplex(dplxr, dna"ATGGCGGGT")  # match 1st barcode without errors
(1, 0)

julia> demultiplex(dplxr, dna"CCGACGGGT")  # match 2nd barcode with one error
(2, 1)
```


[Ref. Sequence demultiplexing](https://biojulia.net/BioSequences.jl/stable/demultiplexer/)


---






# Counting






### Count matches


```julia
julia> count(==, dna"ATTGCA", dna"ATGGCC")
4

julia> matches(dna"ATTGCA", dna"ATGGCC")
4
```


--






### Count mismatches


```julia
julia> count(!=, dna"ATTGCA", dna"ATGGCC")
2

julia> mismatches(dna"ATTGCA", dna"ATGGCC")
2
```


---






### Counting ATCGs


```julia
julia> seq = dna"ATGGCC"
6nt DNA Sequence:
ATGGCC

julia> n = length(seq)
6
```


--


```julia
julia> matches(seq, dna"A"^n)
1

julia> matches(seq, dna"T"^n)
1

julia> matches(seq, dna"C"^n)
2

julia> matches(seq, dna"G"^n)
2
```


--






### Counting ambiguous mers


```julia
julia> count(isambiguous, seq)
0
```


---






# Kmers


  * [KmerAnalysis.jl](https://github.com/BioJulia/KmerAnalysis.jl)
  * [KmerAnalysisMakie.jl](https://github.com/BioJulia/KmerAnalysisMakie.jl)


There are some projects for k-mer analysis but not ready yet...


---






# Alignment


```julia
using BioAlignments
```


--






### General string alignment


```julia
julia> costmodel = CostModel(match=1, mismatch=1, insertion=1, deletion=1);

julia> pairalign(EditDistance(), "book", "cook", costmodel)
PairwiseAlignmentResult{Int64,String,String}:
  distance: 4
  seq: 1 book 4
          |||
  ref: 1 cook 4
```


--






### Distance types


  * `EditDistance`
  * `LevenshteinDistance`
  * `HammingDistance`


---






### Biological sequence alignment


```julia
julia> res = pairalign(GlobalAlignment(), s1, s2, affinegap)
PairwiseAlignmentResult{Int64,LongSequence{DNAAlphabet{4}},LongSequence{DNAAlphabet{4}}}:
  score: -1
  seq:  0 -CCTAG---GA--GGG 10
           ||| |   ||  | |
  ref:  1 ACCTGGTATGATAGCG 16
```


--


```julia
julia> res = pairalign(LocalAlignment(), s1, s2, affinegap)
PairwiseAlignmentResult{Int64,LongSequence{DNAAlphabet{4}},LongSequence{DNAAlphabet{4}}}:
  score: 18
  seq:  1 CCTAGGAGG  9
          ||| | | |
  ref:  2 CCTGGTATG 10
```


--






### Alignment types


  * `GlobalAlignment`: global-to-global alignment
  * `SemiGlobalAlignment`: local-to-global alignment
  * `LocalAlignment`: local-to-local alignment
  * `OverlapAlignment`: end-free alignment


https://biojulia.net/BioAlignments.jl/v0.1/pairalign/


---






# Alignment results


```julia
julia> score(res)
18

julia> aln = alignment(res)
PairwiseAlignment{LongSequence{DNAAlphabet{4}},LongSequence{DNAAlphabet{4}}}:
  seq:  1 CCTAGGAGG  9
          ||| | | |
  ref:  2 CCTGGTATG 10
```


--


```julia
julia> count_matches(aln)
6

julia> count_mismatches(aln)
3
```


--


```julia
julia> count_insertions(aln)
0

julia> count_deletions(aln)
0

julia> count_aligned(aln)
9
```


---






# Get aligned sequences


```julia
julia> collect(aln)
9-element Array{Tuple{DNA,DNA},1}:
 (DNA_C, DNA_C)
 (DNA_C, DNA_C)
 (DNA_T, DNA_T)
 (DNA_A, DNA_G)
 (DNA_G, DNA_G)
 (DNA_G, DNA_T)
 (DNA_A, DNA_A)
 (DNA_G, DNA_T)
 (DNA_G, DNA_G)
```


--


```julia
julia> LongDNASeq([x for (x, _) in aln])
9nt DNA Sequence:
CCTAGGAGG

julia> LongDNASeq([y for (_, y) in aln])
9nt DNA Sequence:
CCTGGTATG
```


---






# Substitution matrices






### For DNA/RNA


```julia
julia> EDNAFULL
SubstitutionMatrix{DNA,Int64}:
     A  C  M  G  R  S  V  T  W  Y  H  K  D  B  N
  A  5 -4  1 -4  1 -4 -1 -4  1 -4 -1 -4 -1 -4 -2
  C -4  5  1 -4 -4  1 -1 -4 -4  1 -1 -4 -4 -1 -2
  M  1  1 -1 -4 -2 -2 -1 -4 -2 -2 -1 -4 -3 -3 -1
  G -4 -4 -4  5  1  1 -1 -4 -4 -4 -4  1 -1 -1 -2
  R  1 -4 -2  1 -1 -2 -1 -4 -2 -4 -3 -2 -1 -3 -1
  S -4  1 -2  1 -2 -1 -1 -4 -4 -2 -3 -2 -3 -1 -1
  V -1 -1 -1 -1 -1 -1 -1 -4 -3 -3 -2 -3 -2 -2 -1
  T -4 -4 -4 -4 -4 -4 -4  5  1  1 -1  1 -1 -1 -2
  W  1 -4 -2 -4 -2 -4 -3  1 -1 -2 -1 -2 -1 -3 -1
  Y -4  1 -2 -4 -4 -2 -3  1 -2 -1 -1 -2 -3 -1 -1
  H -1 -1 -1 -4 -3 -3 -2 -1 -1 -1 -1 -3 -2 -2 -1
  K -4 -4 -4  1 -2 -2 -3  1 -2 -2 -3 -1 -1 -1 -1
  D -1 -4 -3 -1 -1 -3 -2 -1 -1 -3 -2 -1 -1 -2 -1
  B -4 -1 -3 -1 -3 -1 -2 -1 -3 -1 -2 -1 -2 -1 -1
  N -2 -2 -1 -2 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1
(underlined values are default ones)
```


---






### For amino acids


  * PAM: `PAM30`, `PAM70`, `PAM250`
  * BLOSUM: `BLOSUM45`, `BLOSUM50`, `BLOSUM62`, `BLOSUM80`, `BLOSUM90`


--


```julia
julia> BLOSUM62
SubstitutionMatrix{AminoAcid,Int64}:
     A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  O  U  B  J  Z  X  *
  A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0  0̲  0̲ -2  0̲ -1  0 -4
  R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3  0̲  0̲ -1  0̲  0 -1 -4
  N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  0̲  0̲  3  0̲  0 -1 -4
  D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  0̲  0̲  4  0̲  1 -1 -4
  C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1  0̲  0̲ -3  0̲ -3 -2 -4
...
```


---




# I/O






### Read bam file


```julia
using XAM
```


--


```julia
julia> reader = open(BAM.Reader, "data/wgEncodeUwRepliSeqBg02esS1AlnRep1.bam")
XAM.BAM.Reader{IOStream}:
  number of contigs: 24

julia> for record in reader
           if BAM.ismapped(record)
               println(BAM.refname(record), ':', BAM.position(record))
           end
       end
```


--


  * `BAM.Reader`
  * `SAM.Reader`


Data with S1 phase and accession GSM923453 downloaded from [ENCODE](http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeUwRepliSeq).


---






# Getting records in a range


```julia
julia> using GenomicFeatures

julia> reader = open(BAM.Reader, "data/wgEncodeUwRepliSeqBg02esS1AlnRep1.bam", index="data/wgEncodeUwRepliSeqBg02esS1AlnRep1.bam.bai")
XAM.BAM.Reader{IOStream}:
  number of contigs: 24

julia> for record in eachoverlap(reader, "chr1", 53800000:53900000)
           println(BAM.refname(record), ':', BAM.position(record))
       end
chr1:53800902
chr1:53801069
chr1:53801760
chr1:53802091
chr1:53803801
...
```


---






# Getting records that overlaps with genomic features


```julia
using GenomicFeatures, GFF3

features = open(collect, GFF3.Reader, "TAIR10_GFF3_genes.gff")

# Keep mRNA features.
filter!(x -> GFF3.featuretype(x) == "mRNA", features)

# Open a BAM file and iterate over records overlapping mRNA transcripts.
reader = open(BAM.Reader, "SRR1238088.sort.bam", index = "SRR1238088.sort.bam.bai")
for feature in features
    for record in eachoverlap(reader, feature)
        # `record` overlaps `feature`.
        # ...
    end
end
close(reader)
```


Code snipest taken directly from [doc](https://biojulia.net/XAM.jl/stable/man/hts-files/#Getting-records-overlapping-genomic-features-1), but I have no time for that...


---






# GFF file processing


  * [GenomicFeatures.jl documentation](https://biojulia.net/GenomicFeatures.jl/stable/)
  * [GFF3.jl documenttion](https://biojulia.net/GFF3.jl/dev/)


---






# VCF processing


```julia
using GeneticVariation
```


  * [GeneticVariation.jl](https://github.com/BioJulia/GeneticVariation.jl)

      * [Issue](https://github.com/BioJulia/GeneticVariation.jl/issues/32) submitted.
  * [VCFTools.jl](https://github.com/OpenMendel/VCFTools.jl)


---






# Shells and pipelines


  * [Running External Programs - Julia doc](https://docs.julialang.org/en/v1/manual/running-external-programs/)
  * [[Day 5] 在Julia底下用pipeline串連一堆指令](https://ithelp.ithome.com.tw/articles/10200062)


---


class: middle






# Thank you for attention

