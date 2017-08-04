########### Test Suite for Goalign command line tools ##############

set -e
set -u
set -o pipefail


echo "->goalign addid"
cat > expected <<EOF
>prefix_Seq0000_suffix
GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTT
TTTAAACACTTTTACATCGA
>prefix_Seq0001_suffix
TGTCGGACCTAAGTATTGAGTACAACGGTGTATTCCAGCGGTGGAGAGGTCTATTTTTCCGGTTGAAGGACTCTAGAGCT
GTAAAGGGTATGGCCATGTG
>prefix_Seq0002_suffix
CTAAGCGCGGGCGGATTGCTGTTGGAGCAAGGTTAAATACTCGGCAATGCCCCATGATCCCCCAAGGACAATAAGAGCGA
AGTTAGAACAAATGAACCCC
>prefix_Seq0003_suffix
GAGTGGAGGCTTTATGGCACAAGGTATTAGAGACTGAGGGGCACCCCGGCATGGTAAGCAGGAGCCATCGCGAAGGCTTC
AGGTATCTTCCTGTGTTACC
>prefix_Seq0004_suffix
CATAGCCCCTGATGCCCTGACCCGTGTCGCGGCAACGTCTACATTTCACGATAAATACTCCGCTGCTAGTCGGCTCTAGA
TGCTTTTCTTCCAGATCTGG
EOF
goalign random -s 10 -n 5 | goalign addid -n prefix_ | goalign addid -n _suffix -r > result
diff result expected
rm -f expected result mapfile


echo "->goalign clean sites"
cat > expected <<EOF
>Seq0000
ATATGGCGATCAAAGTTCCAATGAGATACTTCCTTTACG
>Seq0001
GCGATTATACTCCGAGGTTTTCCGAATAGTTGTTGGAGT
>Seq0002
TACGTTCGTATAAGCTGCCATCCCGGTGCAGACAATACC
>Seq0003
ATGAGGAAAAATGACGGCTAGCAGCAGGTCGCTCTGTAC
>Seq0004
AACCCCGCCGCACATACGTACTCGCTGTAAGTTCCATTG
>Seq0005
GTGAGCGCTTATGCTTAAGTGATTTGCGTTATATTTCAA
>Seq0006
CCACATTCGTATCTTCCGCCCCCCGGTCCACGCTGCTTC
>Seq0007
TGAACTTTCAGACCAGCGTGCTTTGACACCTGGCTCGGC
>Seq0008
CTAAAACCCTAATCCCTTGCTTCTATCCGCCTGGAGCCT
>Seq0009
CTCGCACGACCAGAGAACTAGGCCGGTTCCGGATATCTA
EOF
cat > expectedlog <<EOF
[Warning] in cmd/clean.go (line 55), message: Alignment (0) length before cleaning=100
[Warning] in cmd/clean.go (line 56), message: Alignment (0) length after cleaning=39
[Warning] in cmd/clean.go (line 57), message: Alignment (0) number of gaps=61
EOF
goalign random -s 10 | goalign mutate gaps -n 1 -r 0.1 -s 10 |  goalign clean sites > result 2>log
diff result expected
rm -f expected result mapfile log expectedlog

echo "->goalign clean seqs"
cat > expected <<EOF
>Seq0000
GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTT
TTTAAACACTTTTACATCGA
>Seq0003
GAGTGGAGGCTTTATGGCACAAGGTATTAGAGACTGAGGGGCACCCCGGCATGGTAAGCAGGAGCCATCGCGAAGGCTTC
AGGTATCTTCCTGTGTTACC
>Seq0005
AGTTTGACTATGAGCGCCGGCTTAGTGCTGACAGTGATGCTCCGTTGTAAGGGTCCTGATGTTCTTGTGCTCGCGCATAT
TAGAGCTGAGTTTCCCAAAG
>Seq0007
CTGGTAATACCTGCGCTATTTCGTCAGTTCGTGTACGGGTAACGATAGCGGTTAATGCTTATTCCGATCAGCTCACACCC
ATGAAGGTGGCTCTGGAGCC
>Seq0009
ACCTACGGCTCTAGACAGCTGAAGTCCGGTTCCGAGCACTGTACGGAAACTTGAAAAGGCTCGACGGAGGCTTGTTCCGC
AGAGTGGGACTATAACATAC
EOF
cat > expectedlog <<EOF
[Warning] in cmd/cleanseqs.go (line 36), message: Alignment (0) #seqs before cleaning=10
[Warning] in cmd/cleanseqs.go (line 37), message: Alignment (0) #seqs after cleaning=5
[Warning] in cmd/cleanseqs.go (line 38), message: Alignment (0) removed sequences=5
EOF
goalign random -s 10 | goalign mutate gaps -n 0.5 -r 0.7 -s 10 |  goalign clean seqs > result 2>log
diff result expected
rm -f expected result mapfile log expectedlog

echo "->goalign random"
cat > expected <<EOF
>Seq0000
GATTAATTTGCCGTAGGCCA
>Seq0001
GAATCTGAAGATCGAACACT
>Seq0002
TTAAGTTTTCACTTCTAATG
>Seq0003
GAGAGGACTAGTTCATACTT
>Seq0004
TTTAAACACTTTTACATCGA
>Seq0005
TGTCGGACCTAAGTATTGAG
>Seq0006
TACAACGGTGTATTCCAGCG
>Seq0007
GTGGAGAGGTCTATTTTTCC
>Seq0008
GGTTGAAGGACTCTAGAGCT
>Seq0009
GTAAAGGGTATGGCCATGTG
EOF
goalign random -s 10 -l 20 > result
diff result expected
rm -f expected result mapfile

echo "->goalign reformat fasta"
cat > expected <<EOF
>Seq0000
GATTA
>Seq0001
ATTTG
>Seq0002
CCGTA
>Seq0003
GGCCA
>Seq0004
GAATC
EOF
goalign random -n 5 -l 5 -s 10 -p | goalign reformat fasta -p > result
diff result expected
rm -f expected result mapfile


echo "->goalign reformat fasta strict"
cat > expected <<EOF
>Seq0000
GATTA
>Seq0001
ATTTG
>Seq0002
CCGTA
>Seq0003
GGCCA
>Seq0004
GAATC
EOF
goalign random -n 5 -l 5 -s 10 -p --input-strict | goalign reformat fasta -p > result
diff result expected
rm -f expected result mapfile


echo "->goalign reformat phylip"
cat > expected <<EOF
   5   5
Seq0000  GATTA
Seq0001  ATTTG
Seq0002  CCGTA
Seq0003  GGCCA
Seq0004  GAATC
EOF
goalign random -n 5 -l 5 -s 10 | goalign reformat phylip > result
diff result expected
rm -f expected result mapfile


echo "->goalign reformat phylip strict"
cat > expected <<EOF
   5   5
Seq0000   GATTA
Seq0001   ATTTG
Seq0002   CCGTA
Seq0003   GGCCA
Seq0004   GAATC
EOF
goalign random -n 5 -l 5 -s 10 | goalign reformat phylip --output-strict > result
diff result expected
rm -f expected result mapfile


echo "->goalign reformat nexus"
cat > expected <<EOF
#NEXUS
begin data;
dimensions ntax=5 nchar=5;
format datatype=nucleotide;
matrix
Seq0000 GATTA
Seq0001 ATTTG
Seq0002 CCGTA
Seq0003 GGCCA
Seq0004 GAATC
;
end;
EOF
goalign random -n 5 -l 5 -s 10 | goalign reformat nexus > result
diff result expected
rm -f expected result mapfile


echo "->goalign reformat tnt"
cat > expected <<EOF
xread

'Tnt input file'

5 5
Seq0000 GATTA
Seq0001 ATTTG
Seq0002 CCGTA
Seq0003 GGCCA
Seq0004 GAATC
;
EOF
goalign random -n 5 -l 5 -s 10 | goalign reformat tnt > result
diff result expected
rm -f expected result mapfile

echo "->goalign compute distance -m f81"
cat > expected <<EOF
5
Tip4	0.000000000000	0.174911845895	0.192803956978	0.232646053483	0.235379041630
Tip0	0.174911845895	0.000000000000	0.082364641962	0.128396525775	0.142776083476
Tip3	0.192803956978	0.082364641962	0.000000000000	0.071285264523	0.086842665158
Tip2	0.232646053483	0.128396525775	0.071285264523	0.000000000000	0.111961817720
Tip1	0.235379041630	0.142776083476	0.086842665158	0.111961817720	0.000000000000
EOF
goalign compute distance -m f81 -i tests/data/test_distance.phy.gz -p > result
diff result expected
rm -f expected result mapfile

echo "->goalign compute distance -m pdist"
cat > expected <<EOF
5
Tip4	0.000000000000	0.156000000000	0.170000000000	0.200000000000	0.202000000000
Tip0	0.156000000000	0.000000000000	0.078000000000	0.118000000000	0.130000000000
Tip3	0.170000000000	0.078000000000	0.000000000000	0.068000000000	0.082000000000
Tip2	0.200000000000	0.118000000000	0.068000000000	0.000000000000	0.104000000000
Tip1	0.202000000000	0.130000000000	0.082000000000	0.104000000000	0.000000000000
EOF
goalign compute distance -m pdist -i tests/data/test_distance.phy.gz -p > result
diff result expected
rm -f expected result mapfile

echo "->goalign compute distance -m jc"
cat > expected <<EOF
5
Tip4	0.000000000000	0.174895415376	0.192783827242	0.232616196228	0.235348439687
Tip0	0.174895415376	0.000000000000	0.082361149505	0.128387859288	0.142765296368
Tip3	0.192783827242	0.082361149505	0.000000000000	0.071282661515	0.086838774745
Tip2	0.232616196228	0.128387859288	0.071282661515	0.000000000000	0.111955277061
Tip1	0.235348439687	0.142765296368	0.086838774745	0.111955277061	0.000000000000
EOF
goalign compute distance -m jc -i tests/data/test_distance.phy.gz -p > result
diff result expected
rm -f expected result mapfile


echo "->goalign compute distance -m k2p"
cat > expected <<EOF
5
Tip4	0.000000000000	0.175065910654	0.193245190192	0.233028942469	0.235571330668
Tip0	0.175065910654	0.000000000000	0.082453892766	0.128434558659	0.142789108437
Tip3	0.193245190192	0.082453892766	0.000000000000	0.071292715047	0.086845484497
Tip2	0.233028942469	0.128434558659	0.071292715047	0.000000000000	0.111966480828
Tip1	0.235571330668	0.142789108437	0.086845484497	0.111966480828	0.000000000000
EOF
goalign compute distance -m k2p -i tests/data/test_distance.phy.gz -p > result
diff result expected
rm -f expected result mapfile


echo "->goalign compute distance -m f84"
cat > expected <<EOF
5
Tip4	0.000000000000	0.175078537680	0.193258659612	0.233050782717	0.235595781888
Tip0	0.175078537680	0.000000000000	0.082456191295	0.128441770138	0.142798708889
Tip3	0.193258659612	0.082456191295	0.000000000000	0.071294953075	0.086849003947
Tip2	0.233050782717	0.128441770138	0.071294953075	0.000000000000	0.111972393604
Tip1	0.235595781888	0.142798708889	0.086849003947	0.111972393604	0.000000000000
EOF
goalign compute distance -m f84 -i tests/data/test_distance.phy.gz -p > result
diff result expected
rm -f expected result mapfile


echo "->goalign compute distance -m tn93"
cat > expected <<EOF
5
Tip4	0.000000000000	0.175562323963	0.193263808045	0.233468379444	0.236008401698
Tip0	0.175562323963	0.000000000000	0.082492772223	0.128715806962	0.142896903086
Tip3	0.193263808045	0.082492772223	0.000000000000	0.071294993201	0.086849459942
Tip2	0.233468379444	0.128715806962	0.071294993201	0.000000000000	0.112213623285
Tip1	0.236008401698	0.142896903086	0.086849459942	0.112213623285	0.000000000000
EOF
goalign compute distance -m tn93 -i tests/data/test_distance.phy.gz -p > result
diff result expected
rm -f expected result mapfile

echo "->goalign compute entropy"
cat > expected <<EOF
Alignment	Site	Entropy
0	0	0.000
0	1	0.000
0	2	0.000
0	3	0.000
0	4	0.000
0	5	0.000
0	6	0.000
0	7	0.000
0	8	0.000
0	9	0.000
0	10	0.500
0	11	0.000
0	12	0.000
0	13	0.500
0	14	0.500
0	15	0.000
0	16	0.000
0	17	0.000
0	18	0.000
EOF
goalign compute entropy -i tests/data/test_distance.phy.gz -p > restmp
head -n 20 restmp > result
diff result expected
rm -f expected result restmp


echo "->goalign compute pssm logo"
cat > expected <<EOF
	A	C	G	T
1	0.004	1.953	0.004	0.004
2	1.953	0.004	0.004	0.004
3	1.953	0.004	0.004	0.004
4	0.004	0.004	0.004	1.953
5	0.004	0.004	1.953	0.004
6	0.004	0.004	0.004	1.953
7	0.004	1.953	0.004	0.004
8	0.004	1.953	0.004	0.004
9	0.004	1.953	0.004	0.004
10	0.004	1.953	0.004	0.004
11	0.002	0.997	0.251	0.002
12	1.953	0.004	0.004	0.004
13	1.953	0.004	0.004	0.004
14	0.002	0.002	0.997	0.251
15	0.002	0.251	0.997	0.002
16	1.953	0.004	0.004	0.004
17	0.004	1.953	0.004	0.004
18	0.004	1.953	0.004	0.004
19	0.004	1.953	0.004	0.004
EOF
goalign compute pssm -n 4 -i tests/data/test_distance.phy.gz -p -c 0.01  > restmp
head -n 20 restmp > result
diff result expected
rm -f expected result restmp


echo "->goalign concat"
cat > expected <<EOF
>Seq0000
GATTAATTTGGATTAATTTG
>Seq0001
CCGTAGGCCACCGTAGGCCA
>Seq0002
GAATCTGAAGGAATCTGAAG
>Seq0003
ATCGAACACTATCGAACACT
>Seq0004
TTAAGTTTTCTTAAGTTTTC
>Seq0005
ACTTCTAATGACTTCTAATG
>Seq0006
GAGAGGACTAGAGAGGACTA
>Seq0007
GTTCATACTTGTTCATACTT
>Seq0008
TTTAAACACTTTTAAACACT
>Seq0009
TTTACATCGATTTACATCGA
EOF
goalign random -l 10 -s 10 | goalign concat  <(goalign random -l 10 -s 10 | goalign shuffle seqs) > result
diff result expected
rm -f expected result

echo "->goalign concat"
cat > expected <<EOF
>Seq0000
GATTAATTTGGATTAATTTG
>Seq0001
CCGTAGGCCACCGTAGGCCA
>Seq0002
GAATCTGAAGGAATCTGAAG
>Seq0003
ATCGAACACTATCGAACACT
>Seq0004
TTAAGTTTTCTTAAGTTTTC
>Seq0005
ACTTCTAATGACTTCTAATG
>Seq0006
GAGAGGACTAGAGAGGACTA
>Seq0007
GTTCATACTTGTTCATACTT
>Seq0008
TTTAAACACTTTTAAACACT
>Seq0009
TTTACATCGATTTACATCGA
EOF
goalign random -l 10 -s 10 | goalign concat  <(goalign random -l 10 -s 10 | goalign shuffle seqs) > result
diff result expected
rm -f expected result

echo "->goalign concat (missing sequences in 1)"
cat > input1 <<EOF
>Seq0000
GACATGGCAG
>Seq0001
AATAGAATAG
>Seq0002
TCGATGTCTT
>Seq0003
CGCGAGCCTC
>Seq0004
AGAAGCTTTA
EOF
cat > input2 <<EOF
>Seq0000
GATTAATTTG
>Seq0001
CCGTAGGCCA
>Seq0002
GAATCTGAAG
>Seq0003
ATCGAACACT
>Seq0004
TTAAGTTTTC
>Seq0005
ACTTCTAATG
>Seq0006
GAGAGGACTA
>Seq0007
GTTCATACTT
>Seq0008
TTTAAACACT
>Seq0009
TTTACATCGA
EOF
cat > expected <<EOF
>Seq0000
GACATGGCAGGATTAATTTG
>Seq0001
AATAGAATAGCCGTAGGCCA
>Seq0002
TCGATGTCTTGAATCTGAAG
>Seq0003
CGCGAGCCTCATCGAACACT
>Seq0004
AGAAGCTTTATTAAGTTTTC
>Seq0005
----------ACTTCTAATG
>Seq0006
----------GAGAGGACTA
>Seq0007
----------GTTCATACTT
>Seq0008
----------TTTAAACACT
>Seq0009
----------TTTACATCGA
EOF
goalign concat -i none input1 input2 > result
diff result expected
rm -f expected result input1 input2

echo "->goalign concat (missing sequences in 2)"
cat > input2 <<EOF
>Seq0000
GACATGGCAG
>Seq0001
AATAGAATAG
>Seq0002
TCGATGTCTT
>Seq0003
CGCGAGCCTC
>Seq0004
AGAAGCTTTA
EOF
cat > input1 <<EOF
>Seq0000
GATTAATTTG
>Seq0001
CCGTAGGCCA
>Seq0002
GAATCTGAAG
>Seq0003
ATCGAACACT
>Seq0004
TTAAGTTTTC
>Seq0005
ACTTCTAATG
>Seq0006
GAGAGGACTA
>Seq0007
GTTCATACTT
>Seq0008
TTTAAACACT
>Seq0009
TTTACATCGA
EOF
cat > expected <<EOF
>Seq0000
GATTAATTTGGACATGGCAG
>Seq0001
CCGTAGGCCAAATAGAATAG
>Seq0002
GAATCTGAAGTCGATGTCTT
>Seq0003
ATCGAACACTCGCGAGCCTC
>Seq0004
TTAAGTTTTCAGAAGCTTTA
>Seq0005
ACTTCTAATG----------
>Seq0006
GAGAGGACTA----------
>Seq0007
GTTCATACTT----------
>Seq0008
TTTAAACACT----------
>Seq0009
TTTACATCGA----------
EOF
goalign concat -i none input1 input2 > result
diff result expected
rm -f expected result input1 input2


echo "->goalign divide"
cat > expected <<EOF
>Seq0000
GATTAATTTG
>Seq0000
GATTAATTTG
>Seq0000
GATTAATTTG
>Seq0000
GATTAATTTG
>Seq0000
GATTAATTTG
>Seq0000
GATTAATTTG
>Seq0000
GATTAATTTG
>Seq0000
GATTAATTTG
>Seq0000
GATTAATTTG
>Seq0000
GATTAATTTG
EOF
rm -f input
for i in {1..10}
do
    goalign random -n 1 -l 10 -s 10 -p >> input
done
goalign divide -i input -p -o divprefix -f
cat divprefix_* > result
diff result expected
rm -f expected result divprefix* input


echo "->goalign mutate gaps"
cat > expected <<EOF
>Seq0000
GATTAATTTGCCGTAGGCCA
>Seq0001
G-ATCTGAAGA-CGAACACT
>Seq0002
TTAAGTTTT-ACT-CTAATG
>Seq0003
GAGAGGACTAGTTCATACTT
>Seq0004
TTTAAACA-TTTTACA-CGA
>Seq0005
TGTCGGACCTAAGTATTGAG
>Seq0006
TAC-ACG-TGTATTCCAGCG
>Seq0007
GTGGAGAGGTCTATTTTTCC
>Seq0008
GGTTGAAGGACT-TAGAGC-
>Seq0009
GTAAAGGGTATGGCCATGTG
EOF
goalign random -s 10 -l 20 | goalign mutate gaps -s 10 > result
diff result expected
rm -f expected result


echo "->goalign mutate snvs"
cat > expected <<EOF
>Seq0000
GATTAATTTCCCGTAGGCCA
>Seq0001
GAATCTGAATATCGAACTAT
>Seq0002
TTAAGTTTTCACTTCTAATG
>Seq0003
GAGAGGACTAGTTCATAATT
>Seq0004
TTTTAACACTTTTACATCGA
>Seq0005
TGTCGGACCTAAGTTTTGTG
>Seq0006
TGCAACGATGTACTCCAGCG
>Seq0007
GTGGAGAGGTCTATTTTTGC
>Seq0008
GGTTAAAGGACTCTATAGCT
>Seq0009
GAAAAGGGTATGGCCATGTG
EOF
goalign random -s 10 -l 20 | goalign mutate snvs -s 10 > result
diff result expected
rm -f expected result


echo "->goalign rename"
cat > mapfile <<EOF
Seq0000	New0000
Seq0001	New0001
Seq0002	New0002
Seq0003	New0003
Seq0004	New0004
EOF
cat > expected <<EOF
>New0000
GATTAATTTGCCGTAGGCCA
>New0001
GAATCTGAAGATCGAACACT
>New0002
TTAAGTTTTCACTTCTAATG
>New0003
GAGAGGACTAGTTCATACTT
>New0004
TTTAAACACTTTTACATCGA
EOF
goalign random -s 10 -l 20 -n 5 | goalign rename -m mapfile > result
diff result expected
rm -f expected result mapfile


echo "->goalign sample seqs"
cat > expected <<EOF
>Seq0001
CCGTAGGCCA
>Seq0002
GAATCTGAAG
>Seq0008
TTTAAACACT
EOF
goalign random -l 10 -s 10 | goalign sample seqs -n 3 -s 10 > result
diff result expected
rm -f expected result



echo "->goalign sample sites"
cat > expected <<EOF
>Seq0000
TTAAT
>Seq0001
GTAGG
>Seq0002
ATCTG
>Seq0003
CGAAC
>Seq0004
AAGTT
>Seq0005
TTCTA
>Seq0006
GAGGA
>Seq0007
TCATA
>Seq0008
TAAAC
>Seq0009
TACAT
EOF
goalign random -l 10 -s 10 | goalign sample sites -l 5 -s 10 > result
diff result expected
rm -f expected result



echo "->goalign shuffle recomb"
cat > expected <<EOF
>Seq0000
GATTAATTTG
>Seq0001
CCGTAGGCCA
>Seq0002
GAATCTGAAG
>Seq0003
ATCGATGAAG
>Seq0004
TTAAGTTTTC
>Seq0005
AAGAGGAATG
>Seq0006
GAGAGGACTA
>Seq0007
GTTCATACTT
>Seq0008
TTTAAACACT
>Seq0009
TTTACATCGA
EOF
goalign random -l 10 -s 10 | goalign shuffle recomb -l 0.5 -n 0.25 -s 11 > result
diff result expected
rm -f expected result

echo "->goalign shuffle recomb"
cat > expected <<EOF
>Seq0000
GATTAATTTG
>Seq0001
CCGTAGGCCA
>Seq0002
GAATCTGAAG
>Seq0003
ATCGATGAAG
>Seq0004
TTAAGTTTTC
>Seq0005
AAGAGGAATG
>Seq0006
GAGAGGACTA
>Seq0007
GTTCATACTT
>Seq0008
TTTAAACACT
>Seq0009
TTTACATCGA
EOF
goalign random -l 10 -s 10 | goalign shuffle recomb -l 0.5 -n 0.25 -s 11 > result
diff result expected
rm -f expected result


echo "->goalign stats"
cat > expected <<EOF
length	10
nseqs	10
avgalleles	3.6000
variable sites	10
char	nb	freq
A	31	0.310000
C	18	0.180000
G	17	0.170000
T	34	0.340000
EOF
goalign random -l 10 -s 10 | goalign stats > result
diff result expected
rm -f expected result


echo "->goalign subseq"
cat > expected <<EOF
>Seq0000
ATTTGCCGTA
>Seq0001
ACTGCGCCGC
>Seq0002
CCATACTCGT
>Seq0003
GCTGTGGAGC
EOF
goalign random -n 4 -s 10 -l 10000 | goalign subseq -l 10 -s 5 > result
diff result expected
rm -f expected result


echo "->goalign subset"
cat > expected <<EOF
>Seq0001
CCGTAGGCCA
>Seq3999
CGGGGCCGAC
EOF
goalign random -n 4000 -s 10 -l 10 | goalign subset Seq0001 Seq3999 > result
diff result expected
rm -f expected result


echo "->goalign trim name"
cat > expected <<EOF
>S01
GATTA
>S02
ATTTG
>S03
CCGTA
>S04
GGCCA
EOF
cat > expectedmap <<EOF
Seq0002	S03
Seq0003	S04
Seq0000	S01
Seq0001	S02
EOF
goalign random -s 10 -n 4 -l 5 | goalign trim name -n 3 -m mapfile > result
diff result expected
diff <(sort mapfile) <(sort expectedmap)
rm -f expected result expectedmap mapfile


echo "->goalign trim seq"
cat > expected <<EOF
>Seq0000
ATTTG
>Seq0001
GGCCA
>Seq0002
TGAAG
>Seq0003
ACACT
EOF
goalign random -s 10 -n 4 -l 10 | goalign trim seq -n 5 -s > result
diff result expected
rm -f expected result 


echo "->goalign unalign"
cat > expected <<EOF
>Seq0000
GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTT
TTTAAACACTTTTACATCGA
>Seq0001
GTCGGACCTAAGATTGAGTACAACGGTGATTCCAGCGGGGAGAGGTCTTTTTTCCGGTGAAGGACTCTAGGCTGTAAAGT
ATGGCATGTG
>Seq0002
CTAACGCGGGCGGATTCTGTGGACAAGGTTAAATACCGGCAATGCCCATGATCCCCCAGGACATAAGAGCGAAGTTAACA
AATGAACCCC
>Seq0003
GAGTGGAGGCTTTATGGCACAAGGTATTAGAGACTGAGGGGCACCCCGGCATGGTAAGCAGGAGCCATCGCGAAGGCTTC
AGGTATCTTCCTGTGTTACC
>Seq0004
CATAGCCCCATGCCGACCGTGTCGCGGCAACGTCTACATTCACGATAAAACTCCGTGCTAGTCGGCTCAGATGCTTTTCT
TCCAATCTGG
>Seq0005
AGTTTGACTATGAGCGCCGGCTTAGTGCTGACAGTGATGCTCCGTTGTAAGGGTCCTGATGTTCTTGTGCTCGCGCATAT
TAGAGCTGAGTTTCCCAAAG
>Seq0006
TCGCACGGGTGGAATGTACGTTATGCATAATCAGCGCTTTCCCGACATGCCCCCTCCGTGGCTCCTTCGCCACGGCGGAC
CTGCGGTTCG
>Seq0007
CTGGTAATACCTGCGCTATTTCGTCAGTTCGTGTACGGGTAACGATAGCGGTTAATGCTTATTCCGATCAGCTCACACCC
ATGAAGGTGGCTCTGGAGCC
>Seq0008
TCGTTAACCCACTCTAACCACCTCTGGCGACATCGGTCTCGGCTTGGATACCTTCTCTATTGGACCCCAGGCTCAACTCG
TGAGCTCTCT
>Seq0009
ACCTACGGCTCTAGACAGCTGAAGTCCGGTTCCGAGCACTGTACGGAAACTTGAAAAGGCTCGACGGAGGCTTGTTCCGC
AGAGTGGGACTATAACATAC
EOF
goalign random -s 10 -p | goalign mutate gaps -s 10 -p | goalign unalign -p > result
diff result expected
rm -f expected result 


echo "->goalign reformat nexus 1"
cat > nexus <<EOF
#NEXUS
BEGIN TAXA;
      TaxLabels fish frog snake mouse;
END;

BEGIN CHARACTERS;
      Dimensions NChar=40;
      Format DataType=DNA;
      Matrix
        fish   ACATA GAGGG TACCT CTAAA
        frog   ACATA GAGGG TACCT CTAAC
        snake  ACATA GAGGG TACCT CTAAG
        mouse  ACATA GAGGG TACCT CTAAT

        fish   ACATA GAGGG TACCT CTAAG
        frog   CCATA GAGGG TACCT CTAAG
        snake  GCATA GAGGG TACCT CTAAG
        mouse  TCATA GAGGG TACCT CTAAG
;
END;

BEGIN TREES;
      Tree best=(fish, (frog, (snake, mouse)));
END;
EOF
cat > expected <<EOF
>fish
ACATAGAGGGTACCTCTAAAACATAGAGGGTACCTCTAAG
>frog
ACATAGAGGGTACCTCTAACCCATAGAGGGTACCTCTAAG
>snake
ACATAGAGGGTACCTCTAAGGCATAGAGGGTACCTCTAAG
>mouse
ACATAGAGGGTACCTCTAATTCATAGAGGGTACCTCTAAG
EOF
goalign reformat fasta -i nexus -x -o result
diff expected result
rm -f expected result nexus


echo "->goalign reformat nexus 2"
cat > nexus <<EOF
#NEXUS
BEGIN TAXA;
      TaxLabels fish frog snake mouse;
END;

BEGIN CHARACTERS;
      Dimensions NChar=40;
      Format DataType=DNA;
      Matrix
        fish   ACATA GAGGG TACCT CTAAA
        fish   ACATA GAGGG TACCT CTAAG

        frog   ACATA GAGGG TACCT CTAAC
        frog   CCATA GAGGG TACCT CTAAG

        snake  ACATA GAGGG TACCT CTAAG
        snake  GCATA GAGGG TACCT CTAAG

        mouse  ACATA GAGGG TACCT CTAAT
        mouse  TCATA GAGGG TACCT CTAAG
;
END;

BEGIN TREES;
      Tree best=(fish, (frog, (snake, mouse)));
END;
EOF
cat > expected <<EOF
>fish
ACATAGAGGGTACCTCTAAAACATAGAGGGTACCTCTAAG
>frog
ACATAGAGGGTACCTCTAACCCATAGAGGGTACCTCTAAG
>snake
ACATAGAGGGTACCTCTAAGGCATAGAGGGTACCTCTAAG
>mouse
ACATAGAGGGTACCTCTAATTCATAGAGGGTACCTCTAAG
EOF
goalign reformat fasta -i nexus -x -o result
diff expected result
rm -f expected result nexus
