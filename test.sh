########### Test Suite for Goalign command line tools ##############

set -e
set -u
set -o pipefail

TESTDATA="tests/data"

GOALIGN=./goalign

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
${GOALIGN} random --seed 10 -n 5 | ${GOALIGN} addid -n prefix_ | ${GOALIGN} addid -n _suffix -r > result
diff -q -b result expected
rm -f expected result mapfile


echo "->goalign addid unaligned"
cat > input <<EOF
>Seq0000
GATTAATTTGCCGTAGGCCAAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTT
TTTAAACACTTTTACATCGA
>Seq0001
TGTCGGACCTAAGTATTGAGGGTCTATTTTTCCGGTTGAAGGACTCTAGAGCT
GTAAAGGGTATGGCCATGTG
>Seq0002
CTAAGCGCGGGCGGATTGCTGTTGGAGCAAGGTTAAATACTCGGCAATGCCCCATGATCCCCCAAGGACAATAAGAGCGA
>Seq0003
AGGTATCTTCCTGTGTTACC
>Seq0004
CA
EOF
cat > expected <<EOF
>prefix_Seq0000_suffix
GATTAATTTGCCGTAGGCCAAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTTTTTAAAC
ACTTTTACATCGA
>prefix_Seq0001_suffix
TGTCGGACCTAAGTATTGAGGGTCTATTTTTCCGGTTGAAGGACTCTAGAGCTGTAAAGGGTATGGCCATGTG
>prefix_Seq0002_suffix
CTAAGCGCGGGCGGATTGCTGTTGGAGCAAGGTTAAATACTCGGCAATGCCCCATGATCCCCCAAGGACAATAAGAGCGA
>prefix_Seq0003_suffix
AGGTATCTTCCTGTGTTACC
>prefix_Seq0004_suffix
CA
EOF
${GOALIGN} addid -n prefix_ -i input --unaligned | ${GOALIGN} addid -n _suffix -r --unaligned > result
diff -q -b result expected
rm -f expected result mapfile input


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
Alignment (0) length before cleaning=100
Alignment (0) length after cleaning=39
Alignment (0) number of gaps=61
Alignment (0) number of start gaps=1
Alignment (0) number of end gaps=1
EOF
${GOALIGN} random --seed 10 | ${GOALIGN} mutate gaps -n 1 -r 0.1 --seed 10 |  ${GOALIGN} clean sites > result 2>log
diff -q -b result expected
diff -q -b log expectedlog
rm -f expected result mapfile log expectedlog


echo "->goalign clean sites --ends"
cat > input <<EOF
>Seq0000
--GGTCCACTCTTTTGTCTT-TACCTA-G-
>Seq0001
G---CACCGGC-CGTAATGACG-ACCC--T
>Seq0002
-T-G-TTTCCTGC-AACAT-ACC-AAC-C-
>Seq0003
A-ACCACAACAGTCA-GTACTCTT-TG--T
>Seq0004
-----GAAGG-CCAAGGT-TCGCCGCCC--
EOF

cat > expected <<EOF
>Seq0000
GTCCACTCTTTTGTCTT-TACCTA
>Seq0001
-CACCGGC-CGTAATGACG-ACCC
>Seq0002
G-TTTCCTGC-AACAT-ACC-AAC
>Seq0003
CCACAACAGTCA-GTACTCTT-TG
>Seq0004
--GAAGG-CCAAGGT-TCGCCGCC
EOF

cat > expectedlog <<EOF
Alignment (0) length before cleaning=30
Alignment (0) length after cleaning=24
Alignment (0) number of gaps=6
Alignment (0) number of start gaps=3
Alignment (0) number of end gaps=3
EOF
${GOALIGN} clean sites -i input -c 0.5 > result 2>log
diff -q -b result expected
diff -q -b log expectedlog
rm -f expected result log expectedlog

echo "->goalign clean sites --ends --char MAJ"
cat > input <<EOF
>A
N-ANGA-GACC
>B
N-TN-T-TTTC
>C
NCTN-TTT--T
>D
N-ANCCCCCCC
EOF

cat > expected <<EOF
>A
AGA-GAC
>B
T-T-TTT
>C
T-TTT--
>D
ACCCCCC
EOF

cat > expected2 <<EOF
>A
ANGA-GAC
>B
TN-T-TTT
>C
TN-TTT--
>D
ANCCCCCC
EOF

cat > expected3 <<EOF
>A
-AGA-GACC
>B
-T-T-TTTC
>C
CT-TTT--T
>D
-ACCCCCCC
EOF


cat > expectedlog <<EOF
Alignment (0) length before cleaning=11
Alignment (0) length after cleaning=7
Alignment (0) number of maj=4
Alignment (0) number of start maj=2
Alignment (0) number of end maj=1
EOF

cat > expectedlog2 <<EOF
Alignment (0) length before cleaning=11
Alignment (0) length after cleaning=8
Alignment (0) number of maj=3
Alignment (0) number of start maj=2
Alignment (0) number of end maj=1
EOF

cat > expectedlog3 <<EOF
Alignment (0) length before cleaning=11
Alignment (0) length after cleaning=9
Alignment (0) number of maj=2
Alignment (0) number of start maj=1
Alignment (0) number of end maj=0
EOF

${GOALIGN} clean sites -i input -c 0.6 --char MAJ > result 2>log
diff -q -b result expected
diff -q -b log expectedlog
rm -f expected result log expectedlog

${GOALIGN} clean sites -i input -c 0.6 --char MAJ --ends > result 2>log
diff -q -b result expected2
diff -q -b log expectedlog2
rm -f expected2 result log expectedlog2

${GOALIGN} clean sites -i input -c 1 --char MAJ > result 2>log
diff -q -b result expected3
diff -q -b log expectedlog3
rm -f expected3 result log expectedlog3

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
${GOALIGN} random --seed 10 | ${GOALIGN} mutate gaps -n 0.5 -r 0.7 --seed 10 |  ${GOALIGN} clean seqs > result 2>log
diff -q -b result expected
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
${GOALIGN} random --seed 10 -l 20 > result
diff -q -b result expected
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
${GOALIGN} random -n 5 -l 5 --seed 10 -p | ${GOALIGN} reformat fasta -p > result
diff -q -b result expected
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
${GOALIGN} random -n 5 -l 5 --seed 10 -p --input-strict | ${GOALIGN} reformat fasta -p > result
diff -q -b result expected
rm -f expected result mapfile

echo "->goalign reformat fasta names with spaces"
cat > input <<EOF
> Seq0000
GATTA
> Seq0001
ATTTG
> Seq0002
CCGTA
> Seq0003
GGCCA
> Seq0004
GAATC
EOF
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
${GOALIGN} reformat fasta -i input -o result
diff -q -b result expected
rm -f expected input result

echo "->goalign reformat fasta sequences with spaces"
cat > input <<EOF
>    Seq0000
GA    T TA
>    Seq0001
A     TT TG
>    Seq0002
CC G TA
>    Seq0003
GG CC A
>    Seq0004
G AA T C
EOF
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
${GOALIGN} reformat fasta -i input -o result
diff -q -b result expected
rm -f expected input result

echo "->goalign reformat phylip"
cat > expected <<EOF
   5   75
Seq0000  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG
Seq0001  TACTTTTTAA ACACTTTTAC ATCGATGTCG GACCTAAGTA TTGAGTACAA CGGTGTATTC
Seq0002  CTATTTTTCC GGTTGAAGGA CTCTAGAGCT GTAAAGGGTA TGGCCATGTG CTAAGCGCGG
Seq0003  AGCAAGGTTA AATACTCGGC AATGCCCCAT GATCCCCCAA GGACAATAAG AGCGAAGTTA
Seq0004  GAGTGGAGGC TTTATGGCAC AAGGTATTAG AGACTGAGGG GCACCCCGGC ATGGTAAGCA

   GAGAGGACTA GTTCA
   CAGCGGTGGA GAGGT
   GCGGATTGCT GTTGG
   GAACAAATGA ACCCC
   GGAGCCATCG CGAAG
EOF
${GOALIGN} random -n 5 -l 75 --seed 10 | ${GOALIGN} reformat phylip > result
diff -q -b result expected
rm -f expected result mapfile


echo "->goalign reformat phylip --one-line"
cat > expected <<EOF
   5   75
Seq0000  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCA
Seq0001  TACTTTTTAA ACACTTTTAC ATCGATGTCG GACCTAAGTA TTGAGTACAA CGGTGTATTC CAGCGGTGGA GAGGT
Seq0002  CTATTTTTCC GGTTGAAGGA CTCTAGAGCT GTAAAGGGTA TGGCCATGTG CTAAGCGCGG GCGGATTGCT GTTGG
Seq0003  AGCAAGGTTA AATACTCGGC AATGCCCCAT GATCCCCCAA GGACAATAAG AGCGAAGTTA GAACAAATGA ACCCC
Seq0004  GAGTGGAGGC TTTATGGCAC AAGGTATTAG AGACTGAGGG GCACCCCGGC ATGGTAAGCA GGAGCCATCG CGAAG
EOF
${GOALIGN} random -n 5 -l 75 --seed 10 | ${GOALIGN} reformat phylip --one-line > result
diff -q -b result expected
rm -f expected result mapfile

echo "->goalign reformat phylip --no-block --one-line"
cat > expected <<EOF
   5   75
Seq0000  GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCA
Seq0001  TACTTTTTAAACACTTTTACATCGATGTCGGACCTAAGTATTGAGTACAACGGTGTATTCCAGCGGTGGAGAGGT
Seq0002  CTATTTTTCCGGTTGAAGGACTCTAGAGCTGTAAAGGGTATGGCCATGTGCTAAGCGCGGGCGGATTGCTGTTGG
Seq0003  AGCAAGGTTAAATACTCGGCAATGCCCCATGATCCCCCAAGGACAATAAGAGCGAAGTTAGAACAAATGAACCCC
Seq0004  GAGTGGAGGCTTTATGGCACAAGGTATTAGAGACTGAGGGGCACCCCGGCATGGTAAGCAGGAGCCATCGCGAAG
EOF
${GOALIGN} random -n 5 -l 75 --seed 10 | ${GOALIGN} reformat phylip --one-line --no-block > result
diff -q -b result expected
rm -f expected result mapfile

echo "->goalign reformat phylip --no-block"
cat > expected <<EOF
   5   75
Seq0000  GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATG
Seq0001  TACTTTTTAAACACTTTTACATCGATGTCGGACCTAAGTATTGAGTACAACGGTGTATTC
Seq0002  CTATTTTTCCGGTTGAAGGACTCTAGAGCTGTAAAGGGTATGGCCATGTGCTAAGCGCGG
Seq0003  AGCAAGGTTAAATACTCGGCAATGCCCCATGATCCCCCAAGGACAATAAGAGCGAAGTTA
Seq0004  GAGTGGAGGCTTTATGGCACAAGGTATTAGAGACTGAGGGGCACCCCGGCATGGTAAGCA

   GAGAGGACTAGTTCA
   CAGCGGTGGAGAGGT
   GCGGATTGCTGTTGG
   GAACAAATGAACCCC
   GGAGCCATCGCGAAG
EOF
${GOALIGN} random -n 5 -l 75 --seed 10 | ${GOALIGN} reformat phylip --no-block > result
diff -q -b result expected
rm -f expected result mapfile


echo "->goalign reformat phylip spaces tabs"
cat > expected <<EOF
   5   5
S-e-q-0-0-00  GATTA
Se-q00-0-1  ATTTG
Se-q00-0-2  CCGTA
Seq00-03  GGCCA
Seq-0004  GAATC
EOF
cat > input.fa <<EOF
> S e q 0	0	00
GATTA
> 	 Se    q00 0 1    	 
ATTTG
>	 Se    q00 0 2    	 	
CCGTA
>		Seq00     03
GGCCA
>Seq		0004
GAATC
EOF
${GOALIGN} reformat phylip -i input.fa --clean-names > result
diff -q -b result expected
rm -f expected result mapfile input.fa


echo "->goalign reformat fasta spaces tabs"
cat > expected <<EOF
>S-e-q-0-0-00-
GATTA
>Se-q00-0-2-
CCGTA
>Seq00-03-
GGCCA
>Se-q-0004
GAATC
EOF
cat > input.fa <<EOF
> S e q 0	0	00[]();.,
GATTA
>	 Se    q00 0 2    	 	[]();.,
CCGTA
>		Seq00     03[]();.,
GGCCA
>Se		q[]();.,0004
GAATC
EOF
${GOALIGN} reformat fasta -i input.fa --clean-names > result
diff -q -b result expected
rm -f expected result mapfile input.fa

echo "->goalign reformat fasta spaces tabs"
cat > expected <<EOF
>seq1-A-simple-comment-to-be-removed?-
GATTA
>seq2-A-second-comment-Maybe-a-wrong-sequence-to-be-updated?-
CCGTA
EOF
cat > input.fa <<EOF
> seq1 [A simple comment; to be removed?]
GATTA
> seq2 [A second comment] (Maybe a wrong sequence: to be updated?)
CCGTA
EOF
${GOALIGN} reformat fasta -i input.fa --clean-names > result
diff -q -b result expected
rm -f expected result mapfile input.fa

echo "->goalign reformat phylip strict"
cat > expected <<EOF
   5   5
Seq0000   GATTA
Seq0001   ATTTG
Seq0002   CCGTA
Seq0003   GGCCA
Seq0004   GAATC
EOF
${GOALIGN} random -n 5 -l 5 --seed 10 | ${GOALIGN} reformat phylip --output-strict > result
diff -q -b result expected
rm -f expected result mapfile

echo "->goalign reformat nexus"
cat > expected <<EOF
#NEXUS
begin data;
dimensions ntax=5 nchar=5;
format datatype=dna;
matrix
Seq0000 GATTA
Seq0001 ATTTG
Seq0002 CCGTA
Seq0003 GGCCA
Seq0004 GAATC
;
end;
EOF
${GOALIGN} random -n 5 -l 5 --seed 10 | ${GOALIGN} reformat nexus > result
diff -q -b result expected
rm -f expected result mapfile

echo "->goalign reformat fasta from clustal"
## Sample file from t-coffee:
## https://github.com/cbcrg/tcoffee/blob/master/examples/sample_aln2.aln
cat > input <<EOF
CLUSTAL W (1.83) multiple sequence alignment

1cms            --GEVASVPLTNYLDSQYFGKIYLGTPPQEFTVLFDTGSSDFWVPSIYCK
4pep            ----IGDEPLENYLDTEYFGTIGIGTPAQDFTVIFDTGSSNLWVPSVYCS
4ape            -STGSATTTPIDSLDDAYITPVQIGTPAQTLNLDFDTGSSDLWVFSSETT
3app            AASGVATNTPTAN-DEEYITPVTIG--GTTLNLNFDTGSADLWVFSTELP
2apr            AGVGTVPMTDYGN-DIEYYGQVTIGTPGKKFNLDFDTGSSDLWIASTLCT
1cms_1          -----------------YTGSLHWVPVTVQQYWQFTVDS---VTISGVVV
                                 *   :            * ..*      *    

1cms            SNACKNHQRFDPR-KSSTFQNLGKPLSIHYGTG-SMQGILGYDTVTVSNI
4pep            SLACSDHNQFNPD-DSSTFEATSQELSITYGTG-SMTGILGYDTVQVGGI
4ape            ASEVDGQTIYTPSKSTTAKLLSGATWSISYGDGSSSSGDVYTDTVSVGGL
3app            ASQQSGHSVYNP--SATGKELSGYTWSISYGDGSSASGNVFTDSVTVGGV
2apr            NCG-SGQTKYDPN-QSSTYQADGRTWSISYGDGSSASGILAKDNVNLGGL
1cms_1          ACEGGCQAILDTG-TSKLVGPSSDILNIQQAIG-ATQNQYGEFDIDCDNL
                      :    .   :.     .   .*  . * :  .      :  ..:

1cms            VDIQQTVGLSTQEPGDVFTYAEFDGILGMAYPSLASEY---SIPVFDNMM
4pep            SDTNQIFGLSETEPGSFLYYAPFDGILGLAYPSISASG---ATPVFDNLW
4ape            TVTGQAVESAKKVSSSFTEDSTIDGLLGLAFSTLNTVSPTQQKTFFDNAK
3app            TAHGQAVQAAQQISAQFQQDTNNDGLLGLAFSSINTVQPQSQTTFFDTVK
2apr            LIKGQTIELAKREAASFASGPN-DGLLGLGFDTITTVR--GVKTPMDNLI
1cms_1          SYMPTVVFEINGKMYPLTPSAYTSQDQGFCTSGFQSEN------HSQKWI
                      .         .   .  .   *:    : :          :.  

1cms            NRHLVAQDLFSVYMDRNGQE--SMLTLGAIDPSY
4pep            DQGLVSQDLFSVYLSSNDDSG-SVVLLGGIDSSY
4ape            AS--LDSPVFTADLGYHAPG---TYNFGFIDTTA
3app            SS--LAQPLFAVALKHQQPG---VYDFGFIDSSK
2apr            SQGLISRPIFGVYLGKAKNGGGGEYIFGGYDSTK
1cms_1          LGDVFIREYYSVFDRANNLVG----LAKAI----
                    .    : .                      


EOF
cat > expected <<EOF
>1cms
--GEVASVPLTNYLDSQYFGKIYLGTPPQEFTVLFDTGSSDFWVPSIYCKSNACKNHQRFDPR-KSSTFQNLGKPLSIHY
GTG-SMQGILGYDTVTVSNIVDIQQTVGLSTQEPGDVFTYAEFDGILGMAYPSLASEY---SIPVFDNMMNRHLVAQDLF
SVYMDRNGQE--SMLTLGAIDPSY
>4pep
----IGDEPLENYLDTEYFGTIGIGTPAQDFTVIFDTGSSNLWVPSVYCSSLACSDHNQFNPD-DSSTFEATSQELSITY
GTG-SMTGILGYDTVQVGGISDTNQIFGLSETEPGSFLYYAPFDGILGLAYPSISASG---ATPVFDNLWDQGLVSQDLF
SVYLSSNDDSG-SVVLLGGIDSSY
>4ape
-STGSATTTPIDSLDDAYITPVQIGTPAQTLNLDFDTGSSDLWVFSSETTASEVDGQTIYTPSKSTTAKLLSGATWSISY
GDGSSSSGDVYTDTVSVGGLTVTGQAVESAKKVSSSFTEDSTIDGLLGLAFSTLNTVSPTQQKTFFDNAKAS--LDSPVF
TADLGYHAPG---TYNFGFIDTTA
>3app
AASGVATNTPTAN-DEEYITPVTIG--GTTLNLNFDTGSADLWVFSTELPASQQSGHSVYNP--SATGKELSGYTWSISY
GDGSSASGNVFTDSVTVGGVTAHGQAVQAAQQISAQFQQDTNNDGLLGLAFSSINTVQPQSQTTFFDTVKSS--LAQPLF
AVALKHQQPG---VYDFGFIDSSK
>2apr
AGVGTVPMTDYGN-DIEYYGQVTIGTPGKKFNLDFDTGSSDLWIASTLCTNCG-SGQTKYDPN-QSSTYQADGRTWSISY
GDGSSASGILAKDNVNLGGLLIKGQTIELAKREAASFASGPN-DGLLGLGFDTITTVR--GVKTPMDNLISQGLISRPIF
GVYLGKAKNGGGGEYIFGGYDSTK
>1cms_1
-----------------YTGSLHWVPVTVQQYWQFTVDS---VTISGVVVACEGGCQAILDTG-TSKLVGPSSDILNIQQ
AIG-ATQNQYGEFDIDCDNLSYMPTVVFEINGKMYPLTPSAYTSQDQGFCTSGFQSEN------HSQKWILGDVFIREYY
SVFDRANNLVG----LAKAI----
EOF
${GOALIGN} reformat fasta -i input -u > result
diff -q -b result expected
rm -f expected result input

echo "->goalign reformat fasta from clustal (2)"
cat > input <<EOF
CLUSTALW

1cms     --GEVASVPLTNYLDSQYFGKIYLGTPPQE 30
4pep     ----IGDEPLENYLDTEYFGTIGIGTPAQD 30
4ape     -STGSATTTPIDSLDDAYITPVQIGTPAQT 30
3app     AASGVATNTPTAN-DEEYITPVTIG--GTT 30
2apr     AGVGTVPMTDYGN-DIEYYGQVTIGTPGKK 30
1cms_1   -----------------YTGSLHWVPVTVQ 30
                          *   :        
EOF
cat > expected <<EOF
>1cms
--GEVASVPLTNYLDSQYFGKIYLGTPPQE
>4pep
----IGDEPLENYLDTEYFGTIGIGTPAQD
>4ape
-STGSATTTPIDSLDDAYITPVQIGTPAQT
>3app
AASGVATNTPTAN-DEEYITPVTIG--GTT
>2apr
AGVGTVPMTDYGN-DIEYYGQVTIGTPGKK
>1cms_1
-----------------YTGSLHWVPVTVQ
EOF
${GOALIGN} reformat fasta -i input -u | sed 's/ (GOALIGN version.*//g'> result
diff -q -b result expected
rm -f expected result input


echo "->goalign reformat clustal from fasta"
cat > input <<EOF
>1cms
--GEVASVPLTNYLDSQYFGKIYLGTPPQE
>4pep
----IGDEPLENYLDTEYFGTIGIGTPAQD
>4ape
-STGSATTTPIDSLDDAYITPVQIGTPAQT
>3app
AASGVATNTPTAN-DEEYITPVTIG--GTT
>2apr
AGVGTVPMTDYGN-DIEYYGQVTIGTPGKK
>1cms_1
-----------------YTGSLHWVPVTVQ
EOF
cat > expected <<EOF
CLUSTAL W

1cms     --GEVASVPLTNYLDSQYFGKIYLGTPPQE 30
4pep     ----IGDEPLENYLDTEYFGTIGIGTPAQD 30
4ape     -STGSATTTPIDSLDDAYITPVQIGTPAQT 30
3app     AASGVATNTPTAN-DEEYITPVTIG--GTT 30
2apr     AGVGTVPMTDYGN-DIEYYGQVTIGTPGKK 30
1cms_1   -----------------YTGSLHWVPVTVQ 30
                          *   :        
EOF
${GOALIGN} reformat clustal -i input | sed 's/ (goalign version.*//g'> result
diff -q -b result expected
rm -f expected result input

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
${GOALIGN} random -n 5 -l 5 --seed 10 | ${GOALIGN} reformat tnt > result
diff -q -b result expected
rm -f expected result mapfile


echo "->goalign reformat tnt auto"
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
${GOALIGN} random -n 5 -l 5 --seed 10 | ${GOALIGN} reformat tnt --auto-detect > result
diff -q -b result expected
rm -f expected result mapfile

echo "->goalign reformat auto from nexus"
cat > expected <<EOF
#NEXUS
begin data;
dimensions ntax=5 nchar=5;
format datatype=dna;
matrix
Seq0000 GATTA
Seq0001 ATTTG
Seq0002 CCGTA
Seq0003 GGCCA
Seq0004 GAATC
;
end;
EOF
${GOALIGN} random -n 5 -l 5 --seed 10 -x | ${GOALIGN} reformat nexus --auto-detect > result
diff -q -b result expected
rm -f expected result mapfile

echo "->goalign reformat auto from phylip strict"
cat > expected <<EOF
   5   5
Seq0000   GATTA
Seq0001   ATTTG
Seq0002   CCGTA
Seq0003   GGCCA
Seq0004   GAATC
EOF
${GOALIGN} random -n 5 -l 5 --seed 10 | ${GOALIGN} reformat phylip --output-strict | ${GOALIGN} reformat phylip --output-strict --auto-detect > result
diff -q -b result expected
rm -f expected result mapfile

echo "->goalign reformat auto from phylip"
cat > expected <<EOF
   5   5
Seq0000  GATTA
Seq0001  ATTTG
Seq0002  CCGTA
Seq0003  GGCCA
Seq0004  GAATC
EOF
${GOALIGN} random -n 5 -l 5 --seed 10 -p | ${GOALIGN} reformat phylip --auto-detect > result
diff -q -b result expected
rm -f expected result mapfile

echo "->goalign reformat auto from fasta"
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
${GOALIGN} random -n 5 -l 5 --seed 10 | ${GOALIGN} reformat fasta --auto-detect > result
diff -q -b result expected
rm -f expected result mapfile

echo "->goalign reformat paml from phylip"
cat > input.test <<EOF
   5   5
Seq0000  GATTA
Seq0001  ATTTG
Seq0002  CCGTA
Seq0003  GGCCA
Seq0004  GAATC
EOF
cat > expected <<EOF
  5 5  I
Seq0000
Seq0001
Seq0002
Seq0003
Seq0004

GATTA
ATTTG
CCGTA
GGCCA
GAATC
EOF
${GOALIGN} reformat paml -i input.test -p > output.paml
diff -q -b output.paml expected
rm -f expected output.paml input.test

echo "->goalign compute distance -m f81"
cat > expected <<EOF
5
Tip4	0.000000000000	0.174911845895	0.192803956978	0.232646053483	0.235379041630
Tip0	0.174911845895	0.000000000000	0.082364641962	0.128396525775	0.142776083476
Tip3	0.192803956978	0.082364641962	0.000000000000	0.071285264523	0.086842665158
Tip2	0.232646053483	0.128396525775	0.071285264523	0.000000000000	0.111961817720
Tip1	0.235379041630	0.142776083476	0.086842665158	0.111961817720	0.000000000000
EOF
${GOALIGN} compute distance -m f81 -i ${TESTDATA}/test_distance.phy.gz -p > result
diff -q -b result expected
rm -f expected result mapfile


echo "->goalign compute distance -m rawdist"
cat > expected <<EOF
5
Tip4	0.000000000000	4.000000000000	3.000000000000	3.000000000000	4.000000000000
Tip0	4.000000000000	0.000000000000	1.000000000000	1.000000000000	2.000000000000
Tip3	3.000000000000	1.000000000000	0.000000000000	0.000000000000	1.000000000000
Tip2	3.000000000000	1.000000000000	0.000000000000	0.000000000000	1.000000000000
Tip1	4.000000000000	2.000000000000	1.000000000000	1.000000000000	0.000000000000
EOF
${GOALIGN} compute distance -m rawdist -i ${TESTDATA}/test_rawdistance.phy.gz -p > result
diff -q -b result expected
${GOALIGN} compute distance -m rawdist --gap-mut 2 -i ${TESTDATA}/test_rawdistance.phy.gz -p > result
diff -q -b result expected
${GOALIGN} compute distance -m rawdist --gap-mut 1 -i ${TESTDATA}/test_rawdistance.phy.gz -p > result
diff -q -b result expected
rm -f expected result mapfile



echo "->goalign compute distance -m rawdist2"
cat > expected <<EOF
5
Tip4	0.000000000000	4.000000000000	3.000000000000	3.000000000000	4.000000000000
Tip0	4.000000000000	0.000000000000	0.000000000000	0.000000000000	1.000000000000
Tip3	3.000000000000	0.000000000000	0.000000000000	0.000000000000	1.000000000000
Tip2	3.000000000000	0.000000000000	0.000000000000	0.000000000000	1.000000000000
Tip1	4.000000000000	1.000000000000	1.000000000000	1.000000000000	0.000000000000
EOF
${GOALIGN} compute distance -m rawdist -i ${TESTDATA}/test_rawdistance2.phy.gz -p > result
diff -q -b result expected
rm -f expected result mapfile


echo "->goalign compute distance -m rawdist2 / gap muts"
cat > expected <<EOF
5
Tip4	0.000000000000	4.000000000000	4.000000000000	4.000000000000	5.000000000000
Tip0	4.000000000000	0.000000000000	1.000000000000	1.000000000000	2.000000000000
Tip3	4.000000000000	1.000000000000	0.000000000000	0.000000000000	1.000000000000
Tip2	4.000000000000	1.000000000000	0.000000000000	0.000000000000	1.000000000000
Tip1	5.000000000000	2.000000000000	1.000000000000	1.000000000000	0.000000000000
EOF
${GOALIGN} compute distance -m rawdist --gap-mut 2 -i ${TESTDATA}/test_rawdistance2.phy.gz -p > result
diff -q -b result expected
${GOALIGN} compute distance -m rawdist --gap-mut 1 -i ${TESTDATA}/test_rawdistance2.phy.gz -p > result
diff -q -b result expected
rm -f expected result mapfile


echo "->goalign compute distance -m pdist / gap muts"
cat > expected <<EOF
5
Tip4	0.000000000000  0.148148148148  0.148148148148  0.148148148148  0.185185185185
Tip0	0.148148148148  0.000000000000  0.037037037037	0.037037037037	0.074074074074
Tip3	0.148148148148	0.037037037037	0.000000000000  0.000000000000  0.038461538462
Tip2	0.148148148148	0.037037037037	0.000000000000  0.000000000000  0.038461538462
Tip1	0.185185185185	0.074074074074	0.038461538462  0.038461538462	0.000000000000
EOF
${GOALIGN} compute distance -m pdist --gap-mut 2 -i ${TESTDATA}/test_rawdistance2.phy.gz -p > result
diff -q -b result expected
${GOALIGN} compute distance -m pdist --gap-mut 1 -i ${TESTDATA}/test_rawdistance2.phy.gz -p > result
diff -q -b result expected
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
${GOALIGN} compute distance -m pdist -i ${TESTDATA}/test_distance.phy.gz -p > result
diff -q -b result expected
${GOALIGN} compute distance -m pdist --gap-mut 2 -i ${TESTDATA}/test_distance.phy.gz -p > result
diff -q -b result expected
${GOALIGN} compute distance -m pdist --gap-mut 1 -i ${TESTDATA}/test_distance.phy.gz -p > result
diff -q -b result expected
rm -f expected result mapfile


echo "->goalign compute distance -m raw / gaps start"
cat >input <<EOF
 5 27
Tip4      CAATGTCCCCCAATCACCCGGTCGACG
Tip0      CAATGTCCCCGAAGGACCCGGTGGACG
Tip3      -------CCCXAAGGACCCG-------
Tip2      ---TGTCCCC-AAGGACCCGGTGGACG
Tip1      CAATGTCCCC-AAGGACCCGGTGGCCG
EOF

cat > expected <<EOF
5
Tip4	0.000000000000	4.000000000000	3.000000000000	4.000000000000	5.000000000000
Tip0	4.000000000000	0.000000000000	1.000000000000	1.000000000000	2.000000000000
Tip3	3.000000000000	1.000000000000	0.000000000000	0.000000000000	0.000000000000
Tip2	4.000000000000	1.000000000000	0.000000000000	0.000000000000	1.000000000000
Tip1	5.000000000000	2.000000000000	0.000000000000	1.000000000000	0.000000000000
EOF
${GOALIGN} compute distance -m rawdist --gap-mut 1 -i input -p > result
diff -q -b result expected
rm -f expected result mapfile

echo "->goalign compute distance -m pdist / gaps start"
cat >input <<EOF
 5 27
Tip4      CAATGTCCCCCAATCACCCGGTCGACG
Tip0      CAATGTCCCCGAAGGACCCGGTGGACG
Tip3      -------CCCXAAGGACCCG-------
Tip2      ---TGTCCCC-AAGGACCCGGTGGACG
Tip1      CAATGTCCCC-AAGGACCCGGTGGCCG
EOF

cat > expected <<EOF
5
Tip4	0.000000000000	0.148148148148	0.230769230769	0.166666666667	0.185185185185
Tip0	0.148148148148	0.000000000000	0.076923076923	0.041666666667	0.074074074074
Tip3	0.230769230769	0.076923076923	0.000000000000	0.000000000000	0.000000000000
Tip2	0.166666666667	0.041666666667	0.000000000000	0.000000000000	0.043478260870
Tip1	0.185185185185	0.074074074074	0.000000000000	0.043478260870	0.000000000000
EOF
${GOALIGN} compute distance -m pdist --gap-mut 1 -i input -p > result
diff -q -b result expected
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
${GOALIGN} compute distance -m jc -i ${TESTDATA}/test_distance.phy.gz -p > result
diff -q -b result expected
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
${GOALIGN} compute distance -m k2p -i ${TESTDATA}/test_distance.phy.gz -p > result
diff -q -b result expected
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
${GOALIGN} compute distance -m f84 -i ${TESTDATA}/test_distance.phy.gz -p > result
diff -q -b result expected
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
${GOALIGN} compute distance -m tn93 -i ${TESTDATA}/test_distance.phy.gz -p > result
diff -q -b result expected
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
${GOALIGN} compute entropy -i ${TESTDATA}/test_distance.phy.gz -p > restmp
head -n 20 restmp > result
diff -q -b result expected
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
${GOALIGN} compute pssm -n 4 -i ${TESTDATA}/test_distance.phy.gz -p -c 0.01  > restmp
head -n 20 restmp > result
diff -q -b result expected
rm -f expected result restmp


echo "->goalign concat 1"
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
${GOALIGN} random -l 10 --seed 10 | ${GOALIGN} concat  <(${GOALIGN} random -l 10 --seed 10 | ${GOALIGN} shuffle seqs) > result
diff -q -b result expected
rm -f expected result

echo "->goalign concat 2"
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
${GOALIGN} random -l 10 --seed 10 | ${GOALIGN} concat  <(${GOALIGN} random -l 10 --seed 10 | ${GOALIGN} shuffle seqs) > result
diff -q -b result expected
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
${GOALIGN} concat -i none input1 input2 > result
diff -q -b result expected
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
${GOALIGN} concat -i none input1 input2 > result
diff -q -b result expected
rm -f expected result input1 input2

echo "->goalign concat (missing sequences in both)"
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
>Seq0000
----------GACATGGCAG
>Seq0001
----------AATAGAATAG
>Seq0002
----------TCGATGTCTT
>Seq0003
----------CGCGAGCCTC
EOF
${GOALIGN} concat -i none input1 input2 > result
diff -q -b result expected
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
    ${GOALIGN} random -n 1 -l 10 --seed 10 -p >> input
done
${GOALIGN} divide -i input -p -o divprefix -f
cat divprefix_* > result
diff -q -b result expected
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
${GOALIGN} random --seed 10 -l 20 | ${GOALIGN} mutate gaps --seed 10 > result
diff -q -b result expected
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
${GOALIGN} random --seed 10 -l 20 | ${GOALIGN} mutate snvs --seed 10 > result
diff -q -b result expected
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
${GOALIGN} random --seed 10 -l 20 -n 5 | ${GOALIGN} rename -m mapfile > result
diff -q -b result expected
rm -f expected result mapfile

echo "->goalign rename regexp"
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
${GOALIGN} random --seed 10 -l 20 -n 5 | ${GOALIGN} rename --regexp 'Seq(\d+)' --replace 'New$1' -m mapfile2 > result
diff -q -b result expected
diff -q -b <(sort mapfile) <(sort mapfile2)
rm -f expected result mapfile mapfile2


echo "->goalign rename --clean-names"
cat > input <<EOF
> S e q 0	0	00[]();.,
GATTA
>	 Se    q00 0 2    	 	[]();.,
CCGTA
>		Seq00     03[]();.,
GGCCA
>Se		q[]();.,0004
GAATC
EOF
cat > expected <<EOF
>S-e-q-0-0-00-
GATTA
>Se-q00-0-2-
CCGTA
>Seq00-03-
GGCCA
>Se-q-0004
GAATC
EOF
cat > expectedmap <<EOF
 S e q 0	0	00[]();.,	S-e-q-0-0-00-
	 Se    q00 0 2    	 	[]();.,	Se-q00-0-2-
		Seq00     03[]();.,	Seq00-03-
Se		q[]();.,0004	Se-q-0004
EOF

${GOALIGN} rename --clean-names -i input -o result --map-file outmap
diff -q -b result expected
rm -f input expected result outmap

echo "->goalign rename --clean-names --unaligned"
cat > input <<EOF
> S e q 0	0	00[]();.,
GATTA
>	 Se    q00 0 2    	 	[]();.,
CCGT
>		Seq00     03[]();.,
GGCCACCCGC
>Se		q[]();.,0004
GA
EOF
cat > expected <<EOF
>S-e-q-0-0-00-
GATTA
>Se-q00-0-2-
CCGT
>Seq00-03-
GGCCACCCGC
>Se-q-0004
GA
EOF
cat > expectedmap <<EOF
 S e q 0	0	00[]();.,	S-e-q-0-0-00-
	 Se    q00 0 2    	 	[]();.,	Se-q00-0-2-
		Seq00     03[]();.,	Seq00-03-
Se		q[]();.,0004	Se-q-0004
EOF

${GOALIGN} rename --unaligned --clean-names -i input -o result --map-file outmap
diff -q -b result expected
rm -f input expected result outmap


echo "->goalign sample seqs"
cat > expected <<EOF
>Seq0001
CCGTAGGCCA
>Seq0002
GAATCTGAAG
>Seq0008
TTTAAACACT
EOF
${GOALIGN} random -l 10 --seed 10 | ${GOALIGN} sample seqs -n 3 --seed 10 > result
diff -q -b result expected
rm -f expected result

echo "->goalign sample seqs unaligned"

cat > input <<EOF
>Seq0000
GATTAATTTG
>Seq0001
CCGTAGG
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
TTTAAACAC
>Seq0009
TTTACATCGA
EOF

cat > expected <<EOF
>Seq0001
CCGTAGG
>Seq0002
GAATCTGAAG
>Seq0008
TTTAAACAC
EOF
${GOALIGN} sample seqs --unaligned -i input -n 3 --seed 10 > result
diff -q -b result expected 
rm -f expected result input



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
${GOALIGN} random -l 10 --seed 10 | ${GOALIGN} sample sites -l 5 --seed 10 > result
diff -q -b result expected
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
${GOALIGN} random -l 10 --seed 10 | ${GOALIGN} shuffle recomb -l 0.5 -n 0.25 --seed 11 > result
diff -q -b result expected
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
${GOALIGN} random -l 10 --seed 10 | ${GOALIGN} shuffle recomb -l 0.5 -n 0.25 --seed 11 > result
diff -q -b result expected
rm -f expected result


echo "->goalign shuffle sites"
cat > expected <<EOF
>Seq0000
GATTGGATCCCGGGGTGCGACACTCGTTCGAGTGTAAATTGCTGCTTTGC
>Seq0001
ACTATTAACGGAGAAGACTAGTTCATGCTTGTGAAGCACTCTGCAAACGA
>Seq0002
TGTATGACGGAAGTACTGCGTACAACTGTGTACTCCAGTGAAACGGAGCT
>Seq0003
CTACGTTTGAGGTTCTAGGACTATAGAGGTGTATACGCCATGGTGACGTG
>Seq0004
CTATACTCCCGCGGTGTGCTGTGGGAGCTAGGATAACGTCTTAACATTCC
>Seq0005
CCCACGATCTCCCACGGAGAATGAGATCAAAGTCATTAGATCGTAACCGC
>Seq0006
CAGGATCGTCTTTTGCGCAAAATGTGAAGGTGCTAGATCGACAACTTGGC
>Seq0007
ATGAGAAGTGGAAAGCATAGGGAAGAGCACTGCAAGAGCCCTTGACCACC
>Seq0008
GTTTTAACCTGGTTAACTCTCCCGTTCTACATTTTTGTATTTGGTTTAGG
>Seq0009
AAATGCGCTCCCCGAGTACCCGACTCGAGAAGTAAACTCTGCAGCTCTTG
EOF

cat > rogueexpected <<EOF
Seq0000
Seq0006
Seq0009
Seq0008
Seq0007
EOF

${GOALIGN} random -l 50 --seed 10 | ${GOALIGN} shuffle sites -r 0.5 --seed 10 --rogue 0.5 --rogue-file rogues.txt > result
diff -q -b result expected
diff -q -b rogueexpected rogues.txt
rm -f expected result rogueexpected rogues.txt

echo "->goalign shuffle sites rogues stables"
cat > expected <<EOF
>Seq0000
CATTTATTTGCGGACGTCGTGAAAGTGTAGTTCAAACACATCAAACTGGC
>Seq0001
CCTTATACCGGAGGATACAAATGGAATAATATTTAGATCTTCAGAAAGCC
>Seq0002
AGTCGCAGCTAAGTATGGGTGATCTCTCTGAATACTATCCGAGTCAAGGG
>Seq0003
CTATGTTTCCGCTGTCGGGGCTCGAGGTCTGTAGAAGGTTTGGCGATAGG
>Seq0004
ATGAGCGCGGGGTAAATGCCCTTTAGACGAAGTCAGCGATTCACCGTTTA
>Seq0005
GCGAGGATCCCCCTCGAAAAATATTAAGGATGTTATCACTATTGTTCATC
>Seq0006
GAATTGAAGCTGTTACGCTATACCGTGGTGGGAATAGAGGGTGACACCCT
>Seq0007
ATAGTAAGCAGTAAGGCTCGCGAATGCCTCAGGTACACTGCTGTCTTCGC
>Seq0008
CACAAGCCTTGCGTGGATCAGCCGGAGCACGGCTACTGCCACTGATCTCC
>Seq0009
TTTACTACTCCACGGCTACACGGACCTAGATGCTTGATTACTAGGTCTGG
EOF

cat > rogueexpected <<EOF
Seq0001
Seq0002
Seq0008
Seq0004
Seq0006
EOF

${GOALIGN} random -l 50 --seed 10 | ${GOALIGN} shuffle sites -r 0.5 --seed 10 --rogue 0.5 --rogue-file rogues.txt --stable-rogues  > result
${GOALIGN} random -l 30 --seed 11 | ${GOALIGN} shuffle sites -r 0.5 --seed 10 --rogue 0.5 --rogue-file rogues2.txt --stable-rogues  > /dev/null
diff -q -b result expected
diff -q -b rogueexpected rogues.txt
# Should be the same list of rogues, even if random gen seed is
# diff -q -berent and length is diff -q -berent (initial seq order is the same)
diff -q -b rogues.txt rogues2.txt
rm -f expected result rogueexpected rogues.txt rogues2.txt

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
alphabet	nucleotide
EOF
${GOALIGN} random -l 10 --seed 10 | ${GOALIGN} stats > result
diff -q -b result expected
rm -f expected result


echo "->goalign stats mutations"
cat > input <<EOF
>A
AAACGACGA-GACC-
>B
--AT-TT-T-TTT--
>C
--CTT-TTT--TCC-
EOF

cat > refseq <<EOF
>ref
CCCCCCCCCCCCCCC
EOF

cat > expected <<EOF
A	9
B	8
C	6
EOF

${GOALIGN} stats mutations -i input --ref-sequence refseq > result
diff -q -b result expected
rm -f expected result refseq

cat > input <<EOF
>A
CCCCCCCCCCCCCCC
>B
--AT-TT-T-TTT--
>C
--CTT-TTT--TCC-
EOF

cat > expected <<EOF
A	0
B	8
C	6
EOF

${GOALIGN} stats mutations -i input --ref-sequence A > result
diff -q -b result expected
rm -f input expected result refseq

echo "->goalign stats gaps"
cat > input <<EOF
>A
AAACGACGA-GACC-
>B
--AT-TT-T-TTT--
>C
--CTT-TTT--TCC-
EOF

cat > expected <<EOF
A	2
B	7
C	6
EOF

${GOALIGN} stats gaps -i input > result
diff -q -b result expected
rm -f expected result

cat > expected <<EOF
A	0
B	2
C	2
EOF

${GOALIGN} stats gaps -i input --from-start > result
diff -q -b result expected
rm -f expected result

cat > expected <<EOF
A	1
B	2
C	1
EOF

${GOALIGN} stats gaps -i input --from-end > result
diff -q -b result expected
rm -f expected result

cat > expected <<EOF
A	0
B	3
C	2
EOF

${GOALIGN} stats gaps -i input --unique > result
diff -q -b result expected
rm -f expected result

cat > expected <<EOF
A	2
B	5
C	4
EOF
${GOALIGN} stats gaps -i input --openning > result
diff -q -b result expected
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
${GOALIGN} random -n 4 --seed 10 -l 10000 | ${GOALIGN} subseq -l 10 -s 5 > result
diff -q -b result expected
rm -f expected result


echo "->goalign subseq / refseq"
cat > input <<EOF
>Seq0000
--ACG--AT-GC
>Seq0001
GGACGTTATCGC
EOF
cat > expected <<EOF
>Seq0000
CG--AT
>Seq0001
CGTTAT
EOF
${GOALIGN} subseq -i input -l 4 -s 1 --ref-seq Seq0000 > result
diff -q -b result expected
rm -f input expected result

echo "->goalign subseq / refseq / rev"
cat > input <<EOF
>Seq0000
--ACG--AT-GC
>Seq0001
GGACGTTATCGC
EOF
cat > expected <<EOF
>Seq0000
--A-GC
>Seq0001
GGACGC
EOF
${GOALIGN} subseq -i input -l 4 -s 1 --ref-seq Seq0000 --reverse > result
diff -q -b result expected
rm -f input expected result


echo "->goalign subseq window phylip"
cat > input <<EOF
   10   5
Seq0000  GGAGT
Seq0001  AAGTA
Seq0002  ACGGT
Seq0003  GCACA
Seq0004  AGTAC
Seq0005  GAAGG
Seq0006  CCTTA
Seq0007  AGCAC
Seq0008  GTACG
Seq0009  GAATA
   10   5
Seq0000  AGAGT
Seq0001  AAGTT
Seq0002  ACGGT
Seq0003  ACACT
Seq0004  AGTAT
Seq0005  AAAGT
Seq0006  ACTTT
Seq0007  AGCAT
Seq0008  ATACT
Seq0009  AAATT
EOF
cat > expected.1 <<EOF
   10   4
Seq0000  GGAG
Seq0001  AAGT
Seq0002  ACGG
Seq0003  GCAC
Seq0004  AGTA
Seq0005  GAAG
Seq0006  CCTT
Seq0007  AGCA
Seq0008  GTAC
Seq0009  GAAT
EOF
cat > expected.2 <<EOF
   10   4
Seq0000  GAGT
Seq0001  AGTA
Seq0002  CGGT
Seq0003  CACA
Seq0004  GTAC
Seq0005  AAGG
Seq0006  CTTA
Seq0007  GCAC
Seq0008  TACG
Seq0009  AATA
EOF
cat > expected.3 <<EOF
   10   4
Seq0000  AGAG
Seq0001  AAGT
Seq0002  ACGG
Seq0003  ACAC
Seq0004  AGTA
Seq0005  AAAG
Seq0006  ACTT
Seq0007  AGCA
Seq0008  ATAC
Seq0009  AAAT
EOF
cat > expected.4 <<EOF
   10   4
Seq0000  GAGT
Seq0001  AGTT
Seq0002  CGGT
Seq0003  CACT
Seq0004  GTAT
Seq0005  AAGT
Seq0006  CTTT
Seq0007  GCAT
Seq0008  TACT
Seq0009  AATT
EOF
${GOALIGN} subseq -i input -p -l 4 -s 0 --step 1 -o output.phylip
diff -q -b output.phylip expected.1
diff -q -b output_sub1.phylip expected.2
diff -q -b output_al1.phylip expected.3
diff -q -b output_al1_sub1.phylip expected.4
rm -f input expected.{1,2,3,4} \
   output.phylip output_sub1.phylip \
   output_al1.phylip output_al1_sub1.phylip

echo "->goalign subset"
cat > expected <<EOF
>Seq0001
CCGTAGGCCA
>Seq3999
CGGGGCCGAC
EOF
${GOALIGN} random -n 4000 --seed 10 -l 10 | ${GOALIGN} subset Seq0001 Seq3999 > result
diff -q -b result expected
rm -f expected result

echo "->goalign subset file"
cat > expected <<EOF
>Seq0001
CCGTAGGCCA
>Seq3999
CGGGGCCGAC
EOF
cat > namefile <<EOF
Seq0001,Seq3999
Seq3999,Seq0001
InexistantSeqName
EOF
${GOALIGN} random -n 4000 --seed 10 -l 10 | ${GOALIGN} subset -f namefile > result
diff -q -b result expected
rm -f expected result namefile

echo "->goalign subset index"
cat > expected <<EOF
>Seq0001
CCGTAGGCCA
>Seq3999
CGGGGCCGAC
EOF
${GOALIGN} random -n 4000 --seed 10 -l 10 | ${GOALIGN} subset --indices 1 3999 > result
diff -q -b result expected
rm -f expected result

echo "->goalign subset index file"
cat > expected <<EOF
>Seq0001
CCGTAGGCCA
>Seq3999
CGGGGCCGAC
EOF
cat > indexfile <<EOF
1,3999
1,3999
100000
EOF
${GOALIGN} random -n 4000 --seed 10 -l 10 | ${GOALIGN} subset --indices -f indexfile > result
diff -q -b result expected
rm -f expected result indexfile

echo "->goalign subset regexps"
cat > expected <<EOF
>Seq0001
CCGTAGGCCA
>Seq3999
CGGGGCCGAC
EOF
${GOALIGN} random -n 4000 --seed 10 -l 10 | ${GOALIGN} subset --regexp '.*0001' '.*3999' > result
diff -q -b result expected
rm -f expected result

echo "->goalign subset regexps file"
cat > expected <<EOF
>Seq0001
CCGTAGGCCA
>Seq3999
CGGGGCCGAC
EOF
cat > regfile <<EOF
.*0001
.*3999
EOF
${GOALIGN} random -n 4000 --seed 10 -l 10 | ${GOALIGN} subset --regexp -f regfile > result
diff -q -b result expected
rm -f expected result regfile


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
${GOALIGN} random --seed 10 -n 4 -l 5 | ${GOALIGN} trim name -n 3 -m mapfile > result
diff -q -b result expected
diff -q -b <(sort mapfile) <(sort expectedmap)
rm -f expected result expectedmap mapfile


echo "->goalign trim name auto"
cat > expected <<EOF
>S1
GATTA
>S2
ATTTG
>S3
CCGTA
>S4
GGCCA
EOF
cat > expectedmap <<EOF
Seq0002	S3
Seq0003	S4
Seq0000	S1
Seq0001	S2
EOF
${GOALIGN} random --seed 10 -n 4 -l 5 | ${GOALIGN} trim name -a -m mapfile > result
diff -q -b result expected
diff -q -b <(sort mapfile) <(sort expectedmap)
rm -f expected result expectedmap mapfile


echo "->goalign trim name auto several ali"
cat > input <<EOF
  4 5
Seq0001 GATTA
Seq0002 ATTTG
Seq0003 CCGTA
Seq0004 GGCCA
  4 5
Seq0005 GATTA
Seq0002 ATTTG
Seq0003 CCGTA
Seq0006 GGCCA
  4 5
Seq0001 GATTA
Seq0006 ATTTG
Seq0007 CCGTA
Seq0005 GGCCA
EOF
cat > expected <<EOF
   4   5
S1  GATTA
S2  ATTTG
S3  CCGTA
S4  GGCCA
   4   5
S5  GATTA
S2  ATTTG
S3  CCGTA
S6  GGCCA
   4   5
S1  GATTA
S6  ATTTG
S7  CCGTA
S5  GGCCA
EOF
cat > mapfile <<EOF
Seq0003	S3
Seq0004	S4
Seq0001	S1
Seq0002	S2
Seq0005	S5
Seq0006	S6
Seq0007	S7
EOF
${GOALIGN} trim name -i input -a -m mapfile2 -p > result
diff -q -b result expected
diff -q -b <(sort mapfile) <(sort mapfile2)
rm -f expected result mapfile input mapfile2


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
${GOALIGN} random --seed 10 -n 4 -l 10 | ${GOALIGN} trim seq -n 5 -s > result
diff -q -b result expected
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
${GOALIGN} random --seed 10 -p | ${GOALIGN} mutate gaps --seed 10 -p | ${GOALIGN} unalign -p > result
diff -q -b result expected
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
${GOALIGN} reformat fasta -i nexus -x -o result
diff -q -b expected result
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
${GOALIGN} reformat fasta -i nexus -x -o result
diff -q -b expected result
rm -f expected result nexus


echo "->goalign sort"
cat > expected <<EOF
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
EOF
${GOALIGN} random --seed 10 -l 10 -n 5 | ${GOALIGN} shuffle seqs | ${GOALIGN} sort > result
diff -q -b expected result
rm -f expected result

echo "->goalign sort --unaligned"
cat > input <<EOF
>Seq0004
TTAAG
>Seq0000
GATTAATT
>Seq0002
GAATCTGAAGC
>Seq0003
ATCGAACACT
>Seq0001
CCGTAGGCCACTG
EOF
cat > expected <<EOF
>Seq0000
GATTAATT
>Seq0001
CCGTAGGCCACTG
>Seq0002
GAATCTGAAGC
>Seq0003
ATCGAACACT
>Seq0004
TTAAG
EOF
${GOALIGN} sort -i input --unaligned > result
diff -q -b expected result
rm -f expected result input

echo "->goalign shuffle seqs --unaligned"
cat > input <<EOF
>Seq0004
TTAAG
>Seq0000
GATTAATT
>Seq0002
GAATCTGAAGC
>Seq0003
ATCGAACACT
>Seq0001
CCGTAGGCCACTG
EOF
cat > expected <<EOF
>Seq0000
GATTAATT
>Seq0001
CCGTAGGCCACTG
>Seq0002
GAATCTGAAGC
>Seq0003
ATCGAACACT
>Seq0004
TTAAG
EOF
${GOALIGN} shuffle seqs -i input --unaligned | ${GOALIGN} sort --unaligned  > result
diff -q -b expected result
rm -f expected result input


echo "->goalign translate"
cat > input <<EOF
>allcodons
GCTGCCGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATG
GATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACT
ACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGA
TAG
EOF
cat > expected <<EOF
>allcodons
AAAALLLLLLRRRRRRKKNNMDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVVV***
EOF
${GOALIGN} translate -i input --phase 0 -o result
diff -q -b expected result
rm -f input expected result

echo "->goalign translate unaligned"
cat > input <<EOF
>allcodons
GCTGCCGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATG
GATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACT
ACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGA
TAG
>allcodons2
GCTGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATG
GATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACT
ACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA
EOF
cat > expected <<EOF
>allcodons
AAAALLLLLLRRRRRRKKNNMDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVVV***
>allcodons2
AAALLLLLLRRRRRRKKNNMDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVV**
EOF
${GOALIGN} translate -i input --phase 0 --unaligned -o result
diff -q -b expected result
rm -f input expected result

echo "->goalign translate 3 phases (standard)"
cat > input <<EOF
>test1 
CATGAGTCTCTCTGATAAGGACAAGGCTGCTGTGAAAGCCCTATGG 
EOF
cat > expected <<EOF
>test1 _0
HESL**GQGCCESPM
>test1 _1
MSLSDKDKAAVKALW
>test1 _2
*VSLIRTRLL*KPY
EOF
${GOALIGN} translate -i input --phase -1 -o result
diff -q -b expected result
rm -f input expected result

echo "->goalign translate 3 phases (mito vertebrate)"
cat > input <<EOF
>s1
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGATAGC
>s2
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA
EOF
cat > expected <<EOF
>s1_0
MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIMVVVV*W*
>s1_1
WMTFSVALPHRNSLPHRVAK*LPQRVAEGGITITLS*LS*CNDS
>s1_2
GWLFLLPSPTATVFLIE*R*DYHNGWR*VASLLHYHSCRSVMM
>s2_0
MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIMVVV*W
>s2_1
WMTFSVALPHRNSLPHRVAK*LPQRVAEGGITITLS*LSYN
>s2_2
GWLFLLPSPTATVFLIE*R*DYHNGWR*VASLLHYHSCRMM
EOF
${GOALIGN} translate -i input --phase -1 --unaligned --genetic-code mitov -o result
diff -q -b expected result
rm -f input expected result

echo "->goalign translate 3 phases (mito invertebrate)"
cat > input <<EOF
>s1
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGATAGC
>s2
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA
EOF
cat > expected <<EOF
>s1_0
MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIMVVVV*W*
>s1_1
WMTFSVALPHRNSLPHRVAK*LPQRVAEGGITITLS*LS*CNDS
>s1_2
GWLFLLPSPTATVFLIE*R*DYHNGWR*VASLLHYHSCRSVMM
>s2_0
MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIMVVV*W
>s2_1
WMTFSVALPHRNSLPHRVAK*LPQRVAEGGITITLS*LSYN
>s2_2
GWLFLLPSPTATVFLIE*R*DYHNGWR*VASLLHYHSCRMM
EOF
${GOALIGN} translate -i input --phase -1 --unaligned --genetic-code mitov -o result
diff -q -b expected result
rm -f input expected result

echo "->goalign translate IUPAC"
cat > input <<EOF
>allcodons
GAYAARUAYCAYAAYUAG
EOF
cat > expected <<EOF
>allcodons
DKYHN*
EOF
${GOALIGN} translate -i input --phase 0 -o result
diff -q -b expected result
rm -f input expected result

echo "->goalign translate IUPAC unaligned"
cat > input <<EOF
>allcodons
GAYAARUAYCAYAAYUAG
EOF
cat > expected <<EOF
>allcodons
DKYHN*
EOF
${GOALIGN} translate -i input --phase 0 --unaligned -o result
diff -q -b expected result
rm -f input expected result

echo "->goalign reformat no ignore identical"
cat > input <<EOF
  6   6
1 AAAAAA
1 AAAAAA
1 AAAAAC
2 GGGGGG
3 GGGGGG
4 GGGGGG
EOF

cat > expected <<EOF
   6   6
1  AAAAAA
1_0001  AAAAAA
1_0002  AAAAAC
2  GGGGGG
3  GGGGGG
4  GGGGGG
EOF
${GOALIGN} reformat phylip -i input -o result -p 
diff -q -b expected result
rm -f input expected result

echo "->goalign reformat ignore identical"
cat > input <<EOF
  6   6
1 AAAAAA
1 AAAAAA
1 AAAAAC
2 GGGGGG
3 GGGGGG
4 GGGGGG
EOF

cat > expected <<EOF
   5   6
1  AAAAAA
1_0001  AAAAAC
2  GGGGGG
3  GGGGGG
4  GGGGGG
EOF
${GOALIGN} reformat phylip -i input -o result --ignore-identical -p 
diff -q -b expected result
rm -f input expected result

echo "->goalign reformat ignore identical 2"
cat > input <<EOF
  6   6
1 AAAAAA
1 AAAAAA
1 AAAAAC
2 GGGGGG
3 GGGGGG
3 GGGGGG
EOF

cat > expected <<EOF
   4   6
1  AAAAAA
1_0001  AAAAAC
2  GGGGGG
3  GGGGGG
EOF
${GOALIGN} reformat phylip -i input -o result --ignore-identical -p 
diff -q -b expected result
rm -f input expected result


echo "->goalign reformat ignore identical clustal->fasta"
cat > input <<EOF
CLUSTAL W (1.83) multiple sequence alignment

1cms            --GEVASVPLTNYLDSQYFGKIYLGTPPQEFTVLFDTGSSDFWVPSIYCK
4pep            ----IGDEPLENYLDTEYFGTIGIGTPAQDFTVIFDTGSSNLWVPSVYCS
4ape            -STGSATTTPIDSLDDAYITPVQIGTPAQTLNLDFDTGSSDLWVFSSETT
1cms            --GEVASVPLTNYLDSQYFGKIYLGTPPQEFTVLFDTGSSDFWVPSIYCK
3app            AASGVATNTPTAN-DEEYITPVTIG--GTTLNLNFDTGSADLWVFSTELP
2apr            AGVGTVPMTDYGN-DIEYYGQVTIGTPGKKFNLDFDTGSSDLWIASTLCT
1cms            --GEVASVPLTNYLDSQYFGKIYLGTPPQEFTVLFDTGSSDFWVPSIYCK
1cms_1          -----------------YTGSLHWVPVTVQQYWQFTVDS---VTISGVVV
                                 *   :            * ..*      *    

1cms            --GEVASVPLTNYLDSQYFGKIYLGTPPQEFTVLFDTGSSDFWVPSIYCK
4pep            ----IGDEPLENYLDTEYFGTIGIGTPAQDFTVIFDTGSSNLWVPSVYCS
4ape            -STGSATTTPIDSLDDAYITPVQIGTPAQTLNLDFDTGSSDLWVFSSETT
1cms            --GEVASVPLTNYLDSQYFGKIYLGTPPQEFTVLFDTGSSDFWVPSIYCK
3app            AASGVATNTPTAN-DEEYITPVTIG--GTTLNLNFDTGSADLWVFSTELP
2apr            AGVGTVPMTDYGN-DIEYYGQVTIGTPGKKFNLDFDTGSSDLWIASTLCT
1cms            --GEVASVPLTNYLDSQYFGKIYLGTPPQEFTVLFDTGSSDFWVPSIYCK
1cms_1          -----------------YTGSLHWVPVTVQQYWQFTVDS---VTISGVVV
                                 *   :            * ..*      *    
EOF

cat > expected <<EOF
>1cms
--GEVASVPLTNYLDSQYFGKIYLGTPPQEFTVLFDTGSSDFWVPSIYCK--GEVASVPLTNYLDSQYFGKIYLGTPPQE
FTVLFDTGSSDFWVPSIYCK
>4pep
----IGDEPLENYLDTEYFGTIGIGTPAQDFTVIFDTGSSNLWVPSVYCS----IGDEPLENYLDTEYFGTIGIGTPAQD
FTVIFDTGSSNLWVPSVYCS
>4ape
-STGSATTTPIDSLDDAYITPVQIGTPAQTLNLDFDTGSSDLWVFSSETT-STGSATTTPIDSLDDAYITPVQIGTPAQT
LNLDFDTGSSDLWVFSSETT
>3app
AASGVATNTPTAN-DEEYITPVTIG--GTTLNLNFDTGSADLWVFSTELPAASGVATNTPTAN-DEEYITPVTIG--GTT
LNLNFDTGSADLWVFSTELP
>2apr
AGVGTVPMTDYGN-DIEYYGQVTIGTPGKKFNLDFDTGSSDLWIASTLCTAGVGTVPMTDYGN-DIEYYGQVTIGTPGKK
FNLDFDTGSSDLWIASTLCT
>1cms_1
-----------------YTGSLHWVPVTVQQYWQFTVDS---VTISGVVV-----------------YTGSLHWVPVTVQ
QYWQFTVDS---VTISGVVV
EOF
${GOALIGN} reformat fasta -u -i input -o result --ignore-identical
diff -q -b expected result
rm -f input expected result

echo "->goalign reformat no ignore identical clustal->fasta"
cat > input <<EOF
CLUSTAL W (1.83) multiple sequence alignment

1cms            --GEVASVPLTNYLDSQYFGKIYLGTPPQEFTVLFDTGSSDFWVPSIYCK
4pep            ----IGDEPLENYLDTEYFGTIGIGTPAQDFTVIFDTGSSNLWVPSVYCS
4ape            -STGSATTTPIDSLDDAYITPVQIGTPAQTLNLDFDTGSSDLWVFSSETT
1cms            --GEVASVPLTNYLDSQYFGKIYLGTPPQEFTVLFDTGSSDFWVPSIYCK
3app            AASGVATNTPTAN-DEEYITPVTIG--GTTLNLNFDTGSADLWVFSTELP
2apr            AGVGTVPMTDYGN-DIEYYGQVTIGTPGKKFNLDFDTGSSDLWIASTLCT
1cms            --GEVASVPLTNYLDSQYFGKIYLGTPPQEFTVLFDTGSSDFWVPSIYCK
1cms_1          -----------------YTGSLHWVPVTVQQYWQFTVDS---VTISGVVV
                                 *   :            * ..*      *    
EOF

cat > expected <<EOF
>1cms
--GEVASVPLTNYLDSQYFGKIYLGTPPQEFTVLFDTGSSDFWVPSIYCK
>4pep
----IGDEPLENYLDTEYFGTIGIGTPAQDFTVIFDTGSSNLWVPSVYCS
>4ape
-STGSATTTPIDSLDDAYITPVQIGTPAQTLNLDFDTGSSDLWVFSSETT
>1cms_0001
--GEVASVPLTNYLDSQYFGKIYLGTPPQEFTVLFDTGSSDFWVPSIYCK
>3app
AASGVATNTPTAN-DEEYITPVTIG--GTTLNLNFDTGSADLWVFSTELP
>2apr
AGVGTVPMTDYGN-DIEYYGQVTIGTPGKKFNLDFDTGSSDLWIASTLCT
>1cms_0002
--GEVASVPLTNYLDSQYFGKIYLGTPPQEFTVLFDTGSSDFWVPSIYCK
>1cms_1
-----------------YTGSLHWVPVTVQQYWQFTVDS---VTISGVVV
EOF
${GOALIGN} reformat fasta -u -i input -o result
diff -q -b expected result
rm -f input expected result

echo "->goalign dedup"
cat > input <<EOF
  5   6
1 AAAAAA
2 CCCCCC
3 GGGGGG
4 GGGGGG
4 GGGGGG
EOF
cat > expected <<EOF
   3   6
1  AAAAAA
2  CCCCCC
3  GGGGGG
EOF
${GOALIGN} dedup -i input -o result -p
diff -q -b expected result
rm -f input expected result

echo "->goalign dedup /2"
cat > input <<EOF
   6   6
3  GGGGGG
1  AAAAAA
1  AAAAAC
2  CCCCCC
4  GGGGGG
4  GGGGGG
EOF
cat > expected <<EOF
   4   6
3  GGGGGG
1  AAAAAA
1_0001  AAAAAC
2  CCCCCC
EOF
${GOALIGN} dedup -i input -o result -p
diff -q -b expected result
rm -f input expected result

echo "->goalign dedup --unaligned"
cat > input <<EOF
>1
AAAAAA
>2
CCCCC
>3
GGGG
>4
GGGGGG
>4
GGGGGG
EOF
cat > expected <<EOF
>1
AAAAAA
>2
CCCCC
>3
GGGG
>4
GGGGGG
EOF
${GOALIGN} dedup -i input -o result --unaligned
diff -q -b expected result
rm -f input expected result

echo "->goalign dedup --unaligned /2"
cat > input <<EOF
>3
GGG
>1
AAAAAA
>1
AAAAA
>2
CC
>4
GGGGGG
>4
GGGGGG
EOF
cat > expected <<EOF
>3
GGG
>1
AAAAAA
>1_0001
AAAAA
>2
CC
>4
GGGGGG
EOF
${GOALIGN} dedup -i input -o result --unaligned
diff -q -b expected result
rm -f input expected result


echo "->goalign build seqboot"
cat > expected.1 <<EOF
>Seq0000
ATTAT
>Seq0001
GTTGT
>Seq0002
ATGAT
>Seq0003
ACCAC
EOF
cat > expected.2 <<EOF
>Seq0000
AAGGT
>Seq0001
GGAAT
>Seq0002
AACCT
>Seq0003
AAGGC
EOF
${GOALIGN} random --seed 10 -l 5 -n 4 -o orig.fa
${GOALIGN} build seqboot --seed 10 -i orig.fa -n 2 -o boot
diff -q -b boot0.fa expected.1
diff -q -b boot1.fa expected.2
diff -q -b orig.fa expected.1 > /dev/null || echo "expected.1 ok"
diff -q -b orig.fa expected.2 > /dev/null || echo "expected.2 ok"
if [[ $(ls boot*.fa| wc -l) -ne 2 ]]; then echo "Wrong number of bootstrap alignments"; exit 1; fi
rm -f boot0.fa boot1.fa expected.2 expected.1 orig.fa

echo "->goalign build seqboot gz"
cat > expected.1 <<EOF
>Seq0000
ATTAT
>Seq0001
GTTGT
>Seq0002
ATGAT
>Seq0003
ACCAC
EOF
cat > expected.2 <<EOF
>Seq0000
AAGGT
>Seq0001
GGAAT
>Seq0002
AACCT
>Seq0003
AAGGC
EOF
${GOALIGN} random --seed 10 -l 5 -n 4 -o orig.fa
${GOALIGN} build seqboot --seed 10 -i orig.fa -n 2 -o boot --gz
diff -q -b <(gunzip -c boot0.fa.gz) expected.1
diff -q -b <(gunzip -c boot1.fa.gz) expected.2
if [[ $(ls boot*.fa.gz| wc -l) -ne 2 ]]; then echo "Wrong number of bootstrap alignments"; exit 1; fi
rm -f boot0.fa.gz boot1.fa.gz expected.2 expected.1 orig.fa

echo "->goalign build seqboot tar gz"
cat > expected.1 <<EOF
>Seq0000
ATTAT
>Seq0001
GTTGT
>Seq0002
ATGAT
>Seq0003
ACCAC
EOF
cat > expected.2 <<EOF
>Seq0000
AAGGT
>Seq0001
GGAAT
>Seq0002
AACCT
>Seq0003
AAGGC
EOF
${GOALIGN} random --seed 10 -l 5 -n 4 -o orig.fa
${GOALIGN} build seqboot --seed 10 -i orig.fa -n 2 -o boot --gz --tar
tar -xzf boot.tar.gz
diff -q -b boot0.fa expected.1
diff -q -b boot1.fa expected.2
if [[ $(ls boot*.fa | wc -l) -ne 2 ]]; then echo "Wrong number of bootstrap alignments"; exit 1; fi
rm -f boot0.fa boot1.fa boot.tar.gz  expected.2 expected.1 orig.fa


echo "->goalign build seqboot partition"
cat > partition <<EOF
M1,part1=1-24/3
M2,part2=2-24/3
M3,part3=3-24/3
EOF
cat > input <<EOF
>Seq0000
ACGACGACGACGACGACGACGACG
>Seq0001
ACGACGACGACGACGACGACGACG
>Seq0002
ACGACGACGACGACGACGACGACG
>Seq0003
ACGACGACGACGACGACGACGACG
EOF
cat > expected <<EOF
>Seq0000
AAAAAAAACCCCCCCCGGGGGGGG
>Seq0001
AAAAAAAACCCCCCCCGGGGGGGG
>Seq0002
AAAAAAAACCCCCCCCGGGGGGGG
>Seq0003
AAAAAAAACCCCCCCCGGGGGGGG
EOF
cat > expected_outpartition <<EOF
M1,part1=1-8
M2,part2=9-16
M3,part3=17-24
EOF

${GOALIGN} build seqboot --seed 10 --partition partition --out-partition out_partition -i input -n 4 -o boot
diff -q -b boot0.fa expected
diff -q -b boot1.fa expected
diff -q -b boot2.fa expected
diff -q -b boot3.fa expected
if [[ $(ls boot*.fa| wc -l) -ne 4 ]]; then echo "Wrong number of bootstrap alignments"; exit 1; fi
rm -f boot0.fa boot1.fa boot2.fa boot3.fa expected input partition out_partition expected_outpartition


echo "->goalign codonalign"
cat > input.aa <<EOF
>Seq0000
D*-AVGQNLK
>Seq0001
IE-FKF-LLM
>Seq0002
ERTSSYFLNT
EOF

cat > input.nt <<EOF
>Seq0000
GATTAAGCCGTAGGCCAGAATCTGAAG
>Seq0001
ATCGAATTTAAGTTTCTTCTAATG
>Seq0002
GAGAGGACTAGTTCATACTTTTTAAACACT
EOF

cat > expected <<EOF
>Seq0000
GATTAA---GCCGTAGGCCAGAATCTGAAG
>Seq0001
ATCGAA---TTTAAGTTT---CTTCTAATG
>Seq0002
GAGAGGACTAGTTCATACTTTTTAAACACT
EOF

${GOALIGN} codonalign -i input.aa -f input.nt -o result
diff -q -b expected result
rm -f expected result input.aa input.nt

echo "->goalign identical"
cat > input1 <<EOF
>Seq0000
GATTAA---GCCGTAGGCCAGAATCTGAAG
>Seq0001
ATCGAA---TTTAAGTTT---CTTCTAATG
>Seq0002
GAGAGGACTAGTTCATACTTTTTAAACACT
EOF

cat > input2 <<EOF
>Seq0001
ATCGAA---TTTAAGTTT---CTTCTAATG
>Seq0000
GATTAA---GCCGTAGGCCAGAATCTGAAG
>Seq0002
GAGAGGACTAGTTCATACTTTTTAAACACT
EOF

cat > input3 <<EOF
>Seq0001
GATTAA---GCCGTAGGCCAGAATCTGAAG
>Seq0000
ATCGAA---TTTAAGTTT---CTTCTAATG
>Seq0002
GAGAGGACTAGTTCATACTTTTTAAACACT
EOF

cat > input4 <<EOF
>Seq0000
GATTAA---GCCGTAGGCCAGAATCTGAAG
>Seq0001
ATCGAA---TTTAAGTTT---CTTCTAATG
>Seq0002
GAGAGGACTAGTTCATACTTTTTAAACACT
>Seq0003
ACGACGACACGACGACGACGAGCAGCAGCA
EOF


cat > expected1 <<EOF
true
EOF

cat > expected2 <<EOF
false
EOF

${GOALIGN} identical -i input1 -c input2 > result
diff -q -b expected1 result

${GOALIGN} identical -i input1 -c input3 > result
diff -q -b expected2 result

${GOALIGN} identical -i input1 -c input4 > result
diff -q -b expected2 result

rm -f input1 input2 input3 input4 expected1 expected2 results


echo "->goalign phase"
cat > input <<EOF
>allcodons
GCTGCCGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATG
GATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACT
ACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGA
TAGC
>allcodons2
GCTGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATG
GATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACT
ACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA
EOF
cat > expected <<EOF
>allcodons
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGG
CGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGATAGC
>allcodons2
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGG
CGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA
EOF

cat > expected.aa <<EOF
>allcodons
MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVVV***
>allcodons2
MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVV**
EOF

${GOALIGN} phase -i input --unaligned -o result --aa-output result.aa
diff -q -b expected result
diff -q -b expected.aa result.aa
rm -f input expected result expected.aa result.aa

echo "->goalign phasent"
cat > input <<EOF
>allcodons
GCTGCCGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATG
GATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACT
ACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGA
TAGC
>allcodons2
GCTGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATG
GATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACT
ACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA
EOF
cat > expected <<EOF
>allcodons
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGG
CGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGATAGC
>allcodons2
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGG
CGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA
EOF

cat > expected.aa <<EOF
>allcodons
MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVVV***
>allcodons2
MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVV**
EOF

${GOALIGN} phasent -i input --unaligned -o result --aa-output result.aa
diff -q -b expected result
diff -q -b expected.aa result.aa
rm -f input expected result expected.aa result.aa


echo "->goalign phase reverse"
cat > input <<EOF
>allcodons
GCTATCATTACACTACGACAACTATGATAATGTAATAGTGATGCCACCCTCCGCCACCCG
TTGTGGTAGTCTCTTCGCTACTCGATGAGGAAGACTGTTGCGGTGGGGGAGGGCAACAGA
AAAAGTCATCCATGTTATTCTTTTTCCTTCTCCGTCGGCGACGCAGTAGGAGAAGCAATA
ACGCTGCGGCAGC
>allcodons2
GCTGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATG
GATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACT
ACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA
EOF
cat > expected <<EOF
>allcodons
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGG
CGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGATAGC
>allcodons2
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGG
CGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA
EOF

${GOALIGN} phase -i input --unaligned -o result  --reverse
diff -q -b expected result
rm -f input expected result

echo "->goalign phasent reverse"
cat > input <<EOF
>allcodons
GCTATCATTACACTACGACAACTATGATAATGTAATAGTGATGCCACCCTCCGCCACCCG
TTGTGGTAGTCTCTTCGCTACTCGATGAGGAAGACTGTTGCGGTGGGGGAGGGCAACAGA
AAAAGTCATCCATGTTATTCTTTTTCCTTCTCCGTCGGCGACGCAGTAGGAGAAGCAATA
ACGCTGCGGCAGC
>allcodons2
GCTGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATG
GATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACT
ACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA
EOF
cat > expected <<EOF
>allcodons
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGG
CGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGATAGC
>allcodons2
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGG
CGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA
EOF
${GOALIGN} phasent -i input --unaligned -o result  --reverse
diff -q -b expected result
rm -f input expected result


echo "->goalign phase reverse cutends"
cat > input <<EOF
>allcodons
GCTATCATTACACTACGACAACTATGATAATGTAATAGTGATGCCACCCTCCGCCACCCG
TTGTGGTAGTCTCTTCGCTACTCGATGAGGAAGACTGTTGCGGTGGGGGAGGGCAACAGA
AAAAGTCATCCATGTTATTCTTTTTCCTTCTCCGTCGGCGACGCAGTAGGAGAAGCAATA
ACGCTGCGGCAGC
>allcodons2
GCTGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATG
GATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACT
ACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA
EOF
cat > expected <<EOF
>allcodons
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGG
CGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAA
>allcodons2
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGG
CGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA
EOF

cat > expected.aa <<EOF
>allcodons
MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVVV*
>allcodons2
MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVV**
EOF

${GOALIGN} phase -i input --unaligned -o result --aa-output result.aa --reverse --cut-end
diff -q -b expected result
diff -q -b expected.aa result.aa
rm -f input expected result expected.aa result.aa

echo "->goalign phasent reverse cutends"
cat > input <<EOF
>allcodons
GCTATCATTACACTACGACAACTATGATAATGTAATAGTGATGCCACCCTCCGCCACCCG
TTGTGGTAGTCTCTTCGCTACTCGATGAGGAAGACTGTTGCGGTGGGGGAGGGCAACAGA
AAAAGTCATCCATGTTATTCTTTTTCCTTCTCCGTCGGCGACGCAGTAGGAGAAGCAATA
ACGCTGCGGCAGC
>allcodons2
GCTGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATG
GATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACT
ACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA
EOF
cat > expected <<EOF
>allcodons
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGG
CGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAA
>allcodons2
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGG
CGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAA
EOF

cat > expected.aa <<EOF
>allcodons
MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVVV*
>allcodons2
MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVV*
EOF

${GOALIGN} phasent -i input --unaligned -o result --aa-output result.aa --reverse --cut-end
diff -q -b expected result
diff -q -b expected.aa result.aa
rm -f input expected result expected.aa result.aa


echo "->goalign sw (smith & waterman)"
cat > input <<EOF
>allcodons
GCTGCCGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATG
GATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACT
ACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGA
TAGC
>allcodons2
GCTGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATG
GATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACT
ACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA
EOF
cat > expected <<EOF
   2   186
allcodons  GCCGCAGCGT TATTGCTTCT CCTACTGCGT CGCCGACGGA GAAGGAAAAA GAATAACATG
allcodons2  GCTGCAGCGT TATTGCTTCT CCTACTGCGT CGCCGACGGA GAAGGAAAAA GAATAACATG

   GATGACTTTT TCTGTTGCCC TCCCCCACCG CAACAGTCTT CCTCATCGAG TAGCGAAGAG
   GATGACTTTT TCTGTTGCCC TCCCCCACCG CAACAGTCTT CCTCATCGAG TAGCGAAGAG

   ACTACCACAA CGGGTGGCGG AGGGTGGCAT CACTATTACA TTATCATAGT TGTCGTAGTG
   ACTACCACAA CGGGTGGCGG AGGGTGGCAT CACTATTACA TTATCATAGT TGTCGTA---

   TAATGA
   TAATGA
EOF

${GOALIGN} sw -i input -o result -p
diff -q -b expected result
rm -f input expected result

echo "->goalign sw / 2"
cat > input <<EOF
>nt1
CTGGGGTTTAACCAGCCATGCCAGTGCAGGTTTAAGAACCGATCCGTACTCTGGGTTACTGATGAAGGATGGGCCGTATC
GCCCCCTTGCGACGTTTCCA
>nt2
TATTATCGTATCGTTTGCATAGACCCGTTATGCCAGCAGATACAGCGTCACAAACTTAGGCTGTAGGGCGTTAGCGGCGC
TCCATGTTTAGACTCACGCC
EOF
cat > expected <<EOF
>nt1
GTTT-----AACCAGCCATGCCAGTGCAGGTTTAAGAACCGATCCGT-----ACTCTGGGTTACTGATGAAGGATGGGCC
GTATCGCCCCCTTGCGACGTTTCCA
>nt2
GTTTGCATAGACCCGTTATGCCA--GCAGAT------ACAG---CGTCACAAACTTAGG----CTG--------TAGGGC
GT---------TAGCGGCG-CTCCA
EOF

${GOALIGN} sw -i input -o result 
diff -q -b expected result
rm -f input expected result


echo "->goalign sw aa / 3"
cat > input <<EOF
>aa1
IDYLPEDDSHMFFTYIFMKNFQALGLWAPLDSVAMLHHQRLHSIRNSARKFVNPEDDAIDYCSLCTYEHVLNNIWNGTSR
YQQIWIKVPQETWPKVKRWM
>aa2
KCDVHGRYDTDREDVSEQTDMPHQRFYSVTSYWWMYQALMGTQALRESQMFAMCWVVCNEQDYKHYYYWEGSTYQYEINQ
GRICKNSVKHNTIGIMRNRI
EOF
cat > expected <<EOF
>aa1
CSLCTYEHVLNNIWNGTSRYQ
>aa2
CNEQDYKHYY--YWEG-STYQ
EOF

${GOALIGN} sw -i input -o result 
diff -q -b expected result
rm -f input expected result


echo "->goalign orf"
cat > input <<EOF
>allcodons
GCTGCCGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATG
GATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACT
ACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGA
TAGC
>allcodons2
GCTGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATG
GATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACT
ACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA
EOF
cat > expected <<EOF
>allcodons
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGG
CGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAA
EOF

${GOALIGN} orf -i input -o result
diff -q -b expected result
rm -f input expected result


echo "->goalign orf reverse"
cat > input <<EOF
>s
GCTATCATTACACTACGACAACTATGATAATGTAATAGTGATGCCACCCTCCGCCACCCGTTGTGGTAGTCTCTTCGCTA
CTCGATGAGGAAGACTGTTGCGGTGGGGGAGGGCAACAGAAAAAGTCATCCATGTTATTCTTTTTCCTTCTCCGTCGGCG
ACGCAGTAGGAGAAGCAATAACGCTGCGGCAGC
EOF

cat > expected <<EOF
>s
ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGG
CGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAA
EOF

${GOALIGN} orf -i input -o result --reverse
diff -q -b expected result
rm -f input expected result

echo "->goalign mask / prot"
cat > input <<EOF
   10   20
Seq0000  PHGVHCVSSY RFEKCPNFFC
Seq0001  EACKWDNTCP MKIETHQHQK
Seq0002  GDMMEDSGSI AIDGIGHHKN
Seq0003  YHPTYLHWSA PDGRCKTQSV
Seq0004  DKVTHHQMWV SDGCFMKQPH
Seq0005  GLKQYYAQRK ATNSHKDLAY
Seq0006  GCRMCILEGL KHVGIQMNWN
Seq0007  NSNLHNANYS QGHVVSDVIF
Seq0008  YYKVPHTLPH QDQSCLTEWG
Seq0009  PGTTAYLLGH DYNWFCSEKN
EOF

cat > expected <<EOF
   10   20
Seq0000  XXGVHCVSSY RFEKCPNFFC
Seq0001  XXCKWDNTCP MKIETHQHQK
Seq0002  XXMMEDSGSI AIDGIGHHKN
Seq0003  XXPTYLHWSA PDGRCKTQSV
Seq0004  XXVTHHQMWV SDGCFMKQPH
Seq0005  XXKQYYAQRK ATNSHKDLAY
Seq0006  XXRMCILEGL KHVGIQMNWN
Seq0007  XXNLHNANYS QGHVVSDVIF
Seq0008  XXKVPHTLPH QDQSCLTEWG
Seq0009  XXTTAYLLGH DYNWFCSEKN
EOF

${GOALIGN} mask -i input -o result -s 0 -l 2 -p
diff -q -b expected result

cat > expected <<EOF
   10   20
Seq0000  PHGVHCVSXX XXEKCPNFFC
Seq0001  EACKWDNTXX XXIETHQHQK
Seq0002  GDMMEDSGXX XXDGIGHHKN
Seq0003  YHPTYLHWXX XXGRCKTQSV
Seq0004  DKVTHHQMXX XXGCFMKQPH
Seq0005  GLKQYYAQXX XXNSHKDLAY
Seq0006  GCRMCILEXX XXVGIQMNWN
Seq0007  NSNLHNANXX XXHVVSDVIF
Seq0008  YYKVPHTLXX XXQSCLTEWG
Seq0009  PGTTAYLLXX XXNWFCSEKN
EOF

${GOALIGN} mask -i input -o result -s 8 -l 4 -p
diff -q -b expected result

cat > expected <<EOF
   10   20
Seq0000  PHGVHCVSSY RFEKCPNFXX
Seq0001  EACKWDNTCP MKIETHQHXX
Seq0002  GDMMEDSGSI AIDGIGHHXX
Seq0003  YHPTYLHWSA PDGRCKTQXX
Seq0004  DKVTHHQMWV SDGCFMKQXX
Seq0005  GLKQYYAQRK ATNSHKDLXX
Seq0006  GCRMCILEGL KHVGIQMNXX
Seq0007  NSNLHNANYS QGHVVSDVXX
Seq0008  YYKVPHTLPH QDQSCLTEXX
Seq0009  PGTTAYLLGH DYNWFCSEXX
EOF

${GOALIGN} mask -i input -o result -s 18 -l 2 -p
diff -q -b expected result

cat > expected <<EOF
   10   20
Seq0000  XXXXXXXXXX XXXXXXXXXX
Seq0001  XXXXXXXXXX XXXXXXXXXX
Seq0002  XXXXXXXXXX XXXXXXXXXX
Seq0003  XXXXXXXXXX XXXXXXXXXX
Seq0004  XXXXXXXXXX XXXXXXXXXX
Seq0005  XXXXXXXXXX XXXXXXXXXX
Seq0006  XXXXXXXXXX XXXXXXXXXX
Seq0007  XXXXXXXXXX XXXXXXXXXX
Seq0008  XXXXXXXXXX XXXXXXXXXX
Seq0009  XXXXXXXXXX XXXXXXXXXX
EOF

${GOALIGN} mask -i input -o result -s 0 -l 20 -p
diff -q -b expected result

cat > expected <<EOF
   10   20
Seq0000  XXXXXXXXXX XXXXXXXXXX
Seq0001  XXXXXXXXXX XXXXXXXXXX
Seq0002  XXXXXXXXXX XXXXXXXXXX
Seq0003  XXXXXXXXXX XXXXXXXXXX
Seq0004  XXXXXXXXXX XXXXXXXXXX
Seq0005  XXXXXXXXXX XXXXXXXXXX
Seq0006  XXXXXXXXXX XXXXXXXXXX
Seq0007  XXXXXXXXXX XXXXXXXXXX
Seq0008  XXXXXXXXXX XXXXXXXXXX
Seq0009  XXXXXXXXXX XXXXXXXXXX
EOF

${GOALIGN} mask -i input -o result -s 0 -l 200 -p
diff -q -b expected result

echo "->goalign mask / nucl"
cat > input <<EOF
   10   20
Seq0000  GATTAATTTG CCGTAGGCCA
Seq0001  GAATCTGAAG ATCGAACACT
Seq0002  TTAAGTTTTC ACTTCTAATG
Seq0003  GAGAGGACTA GTTCATACTT
Seq0004  TTTAAACACT TTTACATCGA
Seq0005  TGTCGGACCT AAGTATTGAG
Seq0006  TACAACGGTG TATTCCAGCG
Seq0007  GTGGAGAGGT CTATTTTTCC
Seq0008  GGTTGAAGGA CTCTAGAGCT
Seq0009  GTAAAGGGTA TGGCCATGTG
EOF

cat > expected <<EOF
   10   20
Seq0000  NNTTAATTTG CCGTAGGCCA
Seq0001  NNATCTGAAG ATCGAACACT
Seq0002  NNAAGTTTTC ACTTCTAATG
Seq0003  NNGAGGACTA GTTCATACTT
Seq0004  NNTAAACACT TTTACATCGA
Seq0005  NNTCGGACCT AAGTATTGAG
Seq0006  NNCAACGGTG TATTCCAGCG
Seq0007  NNGGAGAGGT CTATTTTTCC
Seq0008  NNTTGAAGGA CTCTAGAGCT
Seq0009  NNAAAGGGTA TGGCCATGTG
EOF

${GOALIGN} mask -i input -o result -s 0 -l 2 -p
diff -q -b expected result
rm -f input expected result

echo "->goalign mask / nucl --ref-seq"
cat > input <<EOF
   10   20
Seq0000  GATTAATTTG CCGTAGGCCA
Seq0001  GAATCTGAAG ATCGAACACT
Seq0002  TTAAGT--TC ACTTC--ATG
Seq0003  GAGAGGACTA GTTCATACTT
Seq0004  TTTAAACACT TTTACATCGA
Seq0005  TGTCGGACCT AAGTATTGAG
Seq0006  TACAACGGTG TATTCCAGCG
Seq0007  GTGGAGAGGT CTATTTTTCC
Seq0008  GGTTGAAGGA CTCTAGAGCT
Seq0009  GTAAAGGGTA TGGCCATGTG
EOF

cat > expected <<EOF
   10   20
Seq0000  GATTAATTNN NNNNNNNNCA
Seq0001  GAATCTGANN NNNNNNNNCT
Seq0002  TTAAGT--NN NNNNNNNNTG
Seq0003  GAGAGGACNN NNNNNNNNTT
Seq0004  TTTAAACANN NNNNNNNNGA
Seq0005  TGTCGGACNN NNNNNNNNAG
Seq0006  TACAACGGNN NNNNNNNNCG
Seq0007  GTGGAGAGNN NNNNNNNNCC
Seq0008  GGTTGAAGNN NNNNNNNNCT
Seq0009  GTAAAGGGNN NNNNNNNNTG
EOF

${GOALIGN} mask -i input -o result -s 6 -l 8 --ref-seq Seq0002 -p
diff -q -b expected result
rm -f input expected result

echo "->goalign mask --unique"
cat > input <<EOF
>A
ACANGA-TACC
>B
ACTN-T-TTTC
>C
ACTN-TTT--T
>D
C-ANCCCCCCC
EOF

cat > expected << EOF
>A
ACANNN-TNCC
>B
ACTN-T-TNNC
>C
ACTN-TNT--N
>D
N-ANNNNNNCC
EOF

${GOALIGN} mask --unique -i input -o result
diff -q -b expected result
rm -f input expected result

echo "->goalign replace"
cat > input <<EOF
   10   20
Seq0000  GATTAATTTG CCGTAGGCCA
Seq0001  GAATCTGAAG ATCGAACACT
Seq0002  TTAAGTTTTC ACTTCTAATG
Seq0003  GAGAGGACTA GTTCATACTT
Seq0004  TTTAAACACT TTTACATCGA
Seq0005  TGTCGGACCT AAGTATTGAG
Seq0006  TACAACGGTG TATTCCAGCG
Seq0007  GTGGAGAGGT CTATTTTTCC
Seq0008  GGTTGAAGGA CTCTAGAGCT
Seq0009  GTAAAGGGTA TGGCCATGTG
EOF

cat > expected <<EOF
   10   20
Seq0000  GAT--ATTTG CCG--GGCCA
Seq0001  GAATCTGAAG ATCGAACACT
Seq0002  T--AGTTTTC ACTTC--ATG
Seq0003  GAGAGGAC-- GTTCA--CTT
Seq0004  TT--AACACT TT--CATCGA
Seq0005  TGTCGGACC- -AG--TTGAG
Seq0006  --CAACGGTG --TTCCAGCG
Seq0007  GTGGAGAGGT C--TTTTTCC
Seq0008  GGTTGAAGGA CTC--GAGCT
Seq0009  G--AAGGG-- TGGCCATGTG
EOF

${GOALIGN} replace -s TA -n '--' -i input -o result -p
diff -q -b expected result
rm -f input expected result

echo "->goalign replace unalign"
cat > input <<EOF
>Seq0000
GATTAATTTGCCGTAG
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
GTAAAGGGTATGG
EOF

cat > expected <<EOF
>Seq0000
GAT--ATTTGCCG--G
>Seq0001
GAATCTGAAGATCGAACACT
>Seq0002
T--AGTTTTCACTTC--ATG
>Seq0003
GAGAGGAC--GTTCA--CTT
>Seq0004
TT--AACACTTT--CATCGA
>Seq0005
TGTCGGACC--AG--TTGAG
>Seq0006
--CAACGGTG--TTCCAGCG
>Seq0007
GTGGAGAGGTC--TTTTTCC
>Seq0008
GGTTGAAGGACTC--GAGCT
>Seq0009
G--AAGGG--TGG
EOF

${GOALIGN} replace -s TA -n '--' -i input -o result -p --unaligned
diff -q -b expected result
rm -f input expected result


echo "->goalign replace regex"
cat > input <<EOF
   10   20
Seq0000  GATTAATTTGCCGTAGGCCA
Seq0001  GAATCTGAAGATCGAACACT
Seq0002  TTAAGTTTTCACTTCTAATG
Seq0003  GAGAGGACTAGTTCATACTT
Seq0004  TTTAAACACTTTTACATCGA
Seq0005  TGTCGGACCTAAGTATTGAG
Seq0006  TACAACGGTGTATTCCAGCG
Seq0007  GTGGAGAGGTCTATTTTTCC
Seq0008  GGTTGAAGGACTCTAGAGCT
Seq0009  GTAAAGGGTATGGCCATGTG
EOF

cat > expected <<EOF
   10   20
Seq0000  ---TAATTTG CCGTAGGCCA
Seq0001  ---TCT---- --C---CACT
Seq0002  TTAAGTTTTC ACTTCTAATG
Seq0003  ---AG---TA GTTCATACTT
Seq0004  TTTAAACACT TTTACATCGA
Seq0005  TGTCG---CT AAGTATT---
Seq0006  TACAACGGTG TATTCCAGCG
Seq0007  GTG---AGGT CTATTTTTCC
Seq0008  GGTT---G-- -TCTA---CT
Seq0009  GTAAAGGGTA TGGCCATGTG
EOF

${GOALIGN} replace -s 'GA.' -e -n '---' -p -i input -o result
diff -q -b expected result
rm -f input expected result

echo "->goalign diff"
cat > input <<EOF
   10   100
Seq0000  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0001  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0002  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0003  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0004  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0005  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0006  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0007  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0008  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0009  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
EOF

cat > expected <<EOF
   10   100
Seq0000  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0001  .......... .......... .......... .......... .......... .......... .......... .......... .......... ..........
Seq0002  .......... .......... .......... .......... .......... .......... .......... .......... .......... ..........
Seq0003  .......... .......... .......... .......... .......... .......... .......... .......... .......... ..........
Seq0004  .......... .......... .......... .......... .......... .......... .......... .......... .......... ..........
Seq0005  .......... .......... .......... .......... .......... .......... .......... .......... .......... ..........
Seq0006  .......... .......... .......... .......... .......... .......... .......... .......... .......... ..........
Seq0007  .......... .......... .......... .......... .......... .......... .......... .......... .......... ..........
Seq0008  .......... .......... .......... .......... .......... .......... .......... .......... .......... ..........
Seq0009  .......... .......... .......... .......... .......... .......... .......... .......... .......... ..........
EOF

${GOALIGN} diff -i input -p --one-line -o result
diff -q -b expected result
rm -f input expected result


echo "->goalign diff /2"
cat > input <<EOF
   2   100
Seq0000  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0001  TTAATTTGCC GTAGGCCAGA ATCTGAAGAT CGAACACTTT AAGTTTTCAC TTCTAATGGA GAGGACTAGT TCATACTTTT TAAACACTTT TACATCGAGA
EOF

cat > expected <<EOF
   2   100
Seq0000  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0001  TTAATT.GCC GTAGGCCAG. ATC.GAAG.T CGAAC..TT. AAGTT..CA. TTC.AATGGA ...GACTAGT TCAT.CTT.. .AA.C..TT. .AC.TCGA..
EOF

echo "->goalign diff count"
cat > input <<EOF
   10   100
Seq0000  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0001  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0002  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0003  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0004  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0005  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0006  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0007  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0008  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0009  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
EOF

cat > expected <<EOF
	
Seq0001
Seq0002
Seq0003
Seq0004
Seq0005
Seq0006
Seq0007
Seq0008
Seq0009
EOF

${GOALIGN} diff -i input -p --one-line -o result --counts
diff -q -b expected result
rm -f input expected result

echo "->goalign diff count /2"
cat > input <<EOF
   2   100
Seq0000  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0001  TTAATTTGCC GTAGGCCAGA ATCTGAAGAT CGAACACTTT AAGTTTTCAC TTCTAATGGA GAGGACTAGT TCATACTTTT TAAACACTTT TACATCGAGA
EOF

cat > expected <<EOF
	AC	AG	AT	CA	CG	CT	GA	GC	GT	TA	TC	TG
Seq0001	5	5	12	5	3	7	6	4	4	11	6	6
EOF


${GOALIGN} diff -i input -p --one-line -o result --counts
diff -q -b  expected result
rm -f input expected result


echo "->goalign diff count /3"
cat > input <<EOF
   10   100
Seq0000  GATTAATTTG CCGTAGGCCA GAATCTGAAG ATCGAACACT TTAAGTTTTC ACTTCTAATG GAGAGGACTA GTTCATACTT TTTAAACACT TTTACATCGA
Seq0001  TGTCGGACCT AAGTATTGAG TACAACGGTG TATTCCAGCG GTGGAGAGGT CTATTTTTCC GGTTGAAGGA CTCTAGAGCT GTAAAGGGTA TGGCCATGTG
Seq0002  CTAAGCGCGG GCGGATTGCT GTTGGAGCAA GGTTAAATAC TCGGCAATGC CCCATGATCC CCCAAGGACA ATAAGAGCGA AGTTAGAACA AATGAACCCC
Seq0003  GAGTGGAGGC TTTATGGCAC AAGGTATTAG AGACTGAGGG GCACCCCGGC ATGGTAAGCA GGAGCCATCG CGAAGGCTTC AGGTATCTTC CTGTGTTACC
Seq0004  CATAGCCCCT GATGCCCTGA CCCGTGTCGC GGCAACGTCT ACATTTCACG ATAAATACTC CGCTGCTAGT CGGCTCTAGA TGCTTTTCTT CCAGATCTGG
Seq0005  AGTTTGACTA TGAGCGCCGG CTTAGTGCTG ACAGTGATGC TCCGTTGTAA GGGTCCTGAT GTTCTTGTGC TCGCGCATAT TAGAGCTGAG TTTCCCAAAG
Seq0006  TCGCCACGGT GTGGAATGTA CGTTATGGCA GTAATCAGCG GCTTTCACCG ACATGCCCCC TCCGTGGCTC CTTGCGACCA TCGGCGGACC TGCGGTGTCG
Seq0007  CTGGTAATAC CTGCGCTATT TCGTCAGTTC GTGTACGGGT AACGATAGCG GTTAATGCTT ATTCCGATCA GCTCACACCC ATGAAGGTGG CTCTGGAGCC
Seq0008  TCGTTAACCC ACTCTAACCA CCTCCTGTAG CGACATCGGG TGCTCGGCTT GGATACCTTC GTCATATTGG ACCCCAGGTC TCAACCTCGT GAGCTCTCTG
Seq0009  ACCTACGGCT CTAGACAGCT GAAGTCCGGT TCCGAGCACT GTACGGAAAC TTGAAAAGGC TCGACGGAGG CTTGTTCCGC AGAGTGGGAC TATAACATAC
EOF

cat > expected <<EOF
	AC	AG	AT	CA	CG	CT	GA	GC	GT	TA	TC	TG
Seq0001	5	12	5	5	5	6	2	2	8	7	7	10
Seq0002	5	9	7	6	3	2	3	6	3	13	7	9
Seq0003	4	10	8	5	2	8	3	7	2	6	8	16
Seq0004	8	6	11	5	4	6	1	10	4	7	11	9
Seq0005	8	12	8	5	5	4	4	2	6	7	7	7
Seq0006	10	10	5	3	7	3	3	5	6	3	12	10
Seq0007	6	9	8	2	8	4	2	6	5	9	9	5
Seq0008	11	5	9	3	4	4	4	6	4	5	10	9
Seq0009	7	9	5	4	3	4	4	5	3	9	6	11
EOF

${GOALIGN} diff -i input -p --one-line -o result --counts
diff -q -b  expected result
rm -f input expected result

echo "->goalign diff --reverse"

cat > input <<EOF
>Seq0000
GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTT
TTTAAACACTTTTACATCGA
>Seq0001
TG.CGGACCTAA...TTGAGT.CAAC.GT.TATTCCAG.GG.GGAGAGGTCTA.T.TTCC.GTT.A.GG.C.CT.G.GC.
G.A..GGGTA.GGC...GTG
>Seq0002
CTAAGCGCG.G..G.TTG.T.TTGGA.C.AGGTT..ATAC.CGGCAA.G.C.CATG.TCCCCC.A.GAC.A.AAGAG.GA
AG.T.GA..AAA.GA.C.CC
>Seq0003
..G.GGAGGCTTTAT...ACA.GGTATT...GACTGAGGGGC.CCCCGG..TGGTA.GCA.GAGCC.TCGCGAAGGCT.C
AGGT.T.TTCC.GTGT.ACC
>Seq0004
C..AGCCCCTGATGCCCTG.CCCGTGTCGCGG.A.CGT..AC.TT.CACG.TAAA..C.CCGCT.CTAGTCGG.TCTAGA
.GCTTTTCT.CCAGATCT.G
>Seq0005
AG..TGAC.ATGAGC.C.GGCTTAG..CT..CA.TGATGC.CCGT.G.AAGGG..CTGAT.TTCTTGTGCTCG.GC.TA.
.AG.GCTGAG...C.CAAAG
>Seq0006
TCGCC.CGGTGT.G.ATGT.CGT.A..GCAG.AATCAG.GGCTTTCACCG..A.GCCCCCTCCGT.G..CC..GCG..CA
.CGGCGG..C.GCGGTGTCG
>Seq0007
CTGGT.A.AC.T.CGCTATTTCG..A.TTCG.GT.CGGG.AACGA.AGCGGT.AA.GC.TATTCC..TC..C...C..CC
A.G..GGTGGC.CTGGAGCC
>Seq0008
TCG.T.ACCCA.TCTAA...CCTC...T..CGAC.T.GGG.GCTCGGC.TGGA.ACCT.C.TC.TATTGGACC.CAGG.C
.CA.CCTCG.GAGCTC..TG
>Seq0009
ACC..CGGCT.TAG.CAG.T...GTCCGGTTC...G....G..C.GAAA.TTGAAA.GGCTC..C.GAGGC..GT.C.GC
AGAGTGGGAC.A..ACATAC
EOF

cat > expected <<EOF
>Seq0000
GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTT
TTTAAACACTTTTACATCGA
>Seq0001
TGTCGGACCTAAGTATTGAGTACAACGGTGTATTCCAGCGGTGGAGAGGTCTATTTTTCCGGTTGAAGGACTCTAGAGCT
GTAAAGGGTATGGCCATGTG
>Seq0002
CTAAGCGCGGGCGGATTGCTGTTGGAGCAAGGTTAAATACTCGGCAATGCCCCATGATCCCCCAAGGACAATAAGAGCGA
AGTTAGAACAAATGAACCCC
>Seq0003
GAGTGGAGGCTTTATGGCACAAGGTATTAGAGACTGAGGGGCACCCCGGCATGGTAAGCAGGAGCCATCGCGAAGGCTTC
AGGTATCTTCCTGTGTTACC
>Seq0004
CATAGCCCCTGATGCCCTGACCCGTGTCGCGGCAACGTCTACATTTCACGATAAATACTCCGCTGCTAGTCGGCTCTAGA
TGCTTTTCTTCCAGATCTGG
>Seq0005
AGTTTGACTATGAGCGCCGGCTTAGTGCTGACAGTGATGCTCCGTTGTAAGGGTCCTGATGTTCTTGTGCTCGCGCATAT
TAGAGCTGAGTTTCCCAAAG
>Seq0006
TCGCCACGGTGTGGAATGTACGTTATGGCAGTAATCAGCGGCTTTCACCGACATGCCCCCTCCGTGGCTCCTTGCGACCA
TCGGCGGACCTGCGGTGTCG
>Seq0007
CTGGTAATACCTGCGCTATTTCGTCAGTTCGTGTACGGGTAACGATAGCGGTTAATGCTTATTCCGATCAGCTCACACCC
ATGAAGGTGGCTCTGGAGCC
>Seq0008
TCGTTAACCCACTCTAACCACCTCCTGTAGCGACATCGGGTGCTCGGCTTGGATACCTTCGTCATATTGGACCCCAGGTC
TCAACCTCGTGAGCTCTCTG
>Seq0009
ACCTACGGCTCTAGACAGCTGAAGTCCGGTTCCGAGCACTGTACGGAAACTTGAAAAGGCTCGACGGAGGCTTGTTCCGC
AGAGTGGGACTATAACATAC
EOF

${GOALIGN} diff --reverse -i input -o result
diff -q -b  expected result
rm -f input expected result


echo "->goalign stats auto from 1000 random phylip alignments"
cat > expected.nb <<EOF
100
EOF
cat > expected.len <<EOF
5
EOF
for i in {1..100}; do ${GOALIGN} random -p -n 10 -l 5; done | ${GOALIGN} stats nseq --auto-detect | wc -l | sed 's/^[ \t]*//g' > result.nb
for i in {1..100}; do ${GOALIGN} random -p -n 10 -l 5; done | ${GOALIGN} stats length --auto-detect | sort -u > result.len

diff -q -b result.nb expected.nb
diff -q -b result.len expected.len

rm -f expected.len expected.nb result.len result.nb

echo "->goalign nseq"
cat > input <<EOF
>1
GGGGGGGGGGGGGGG
>2
TTTTTTTTTTTTTTTTTTTT
>3
GGGGGG
>4
AAAAAAAAAA
EOF

cat > expected <<EOF
4
EOF
${GOALIGN} stats nseq -i input --unaligned > result
diff -q -b  expected result
rm -f input expected result

echo "->goalign compress"

cat > input <<EOF
   4   20
1 GGGGGGGGGGGGGGGGGGGG
2 TTTTTTTTTTTTTTTTTTTT
3 GGGGGGGGGGCCCCTTTTTT
4 AAAAAAAAAAAAAAAAAAAA
EOF

cat > expected <<EOF
   4   3
1 GGG
2 TTT
3 CGT
4 AAA
EOF

cat > wexp <<EOF
4
10
6
EOF

${GOALIGN} compress -i input -p -o result --weight-out wres
diff -q -b  expected result
diff -q -b  wres wexp

rm -f input expected result wexp wres



echo "->goalign split"
cat > input <<EOF
>s1
AAAACCCCCGG
>2
AAAACCCCCGG
>3
AAAACCCCCGG
>4
AAAACCCCCGG
>5
AAAACCCCCGG
EOF
cat > partitions <<EOF
M1,p1=1-4,10-11
M2,p2=5-9
EOF

cat > exp_p1 <<EOF
>s1
AAAAGG
>2
AAAAGG
>3
AAAAGG
>4
AAAAGG
>5
AAAAGG
EOF

cat > exp_p2 <<EOF
>s1
CCCCC
>2
CCCCC
>3
CCCCC
>4
CCCCC
>5
CCCCC
EOF

${GOALIGN} split -i input --partition partitions --out-prefix ./
diff -q -b exp_p1 p1.fa
diff -q -b exp_p2 p2.fa
rm -f input exp_p1 exp_p2 partitions p1.fa p2.fa

echo "->goalign split codons/2"
cat > input <<EOF
>s1
ACGACGACGACG
>2
ACGACGACGACG
>3
ACGACGACGACG
>4
ACGACGACGACG
>5
CGACGACGACGA
EOF
cat > partitions <<EOF
M1,p1=1-12/3,2-12/3
M2,p2=3-12/3
EOF

cat > exp_p1 <<EOF
>s1
ACACACAC
>2
ACACACAC
>3
ACACACAC
>4
ACACACAC
>5
CGCGCGCG
EOF

cat > exp_p2 <<EOF
>s1
GGGG
>2
GGGG
>3
GGGG
>4
GGGG
>5
AAAA
EOF

${GOALIGN} split -i input --partition partitions --out-prefix ./
diff -q -b exp_p1 p1.fa
diff -q -b exp_p2 p2.fa
rm -f input exp_p1 exp_p2 partitions p1.fa p2.fa


echo "->goalign consensus"
cat > input <<EOF
>s1
ACGACGACGACG
>2
ACGACGACGACG
>3
TTTTTTTTTTTT
EOF
cat > expected <<EOF
>consensus
ACGACGACGACG
EOF

${GOALIGN} consensus -i input -o result
diff -q -b expected result
rm -f input expected result

echo "->goalign consensus /2"
cat > input <<EOF
>s1
ACGACGACGACC
>2
ATCTT-TTTTTC
>3
ATCTT-TTTTTT
EOF
cat > expected <<EOF
>consensus
ATCTT-TTTTTC
EOF

${GOALIGN} consensus -i input -o result
diff -q -b expected result
rm -f input expected result


echo "->goalign append"
cat > input.1 <<EOF
>s1
ACGACGACGACC
>2
ATCTT-TTTTTC
>3
ATCTT-TTTTTT
EOF

cat > input.2 <<EOF
>s4
ACGACGACGACC
>5
ATCTT-TTTTTC
>6
ATCTT-TTTTTT
EOF

cat > input.3 <<EOF
>s7
ACGACGACGACC
>8
ATCTT-TTTTTC
>9
ATCTT-TTTTTT
EOF

cat > expected <<EOF
>s1
ACGACGACGACC
>2
ATCTT-TTTTTC
>3
ATCTT-TTTTTT
>s4
ACGACGACGACC
>5
ATCTT-TTTTTC
>6
ATCTT-TTTTTT
>s7
ACGACGACGACC
>8
ATCTT-TTTTTC
>9
ATCTT-TTTTTT
EOF

${GOALIGN} append -i input.1 input.2 input.3 -o result
diff -q -b expected result
rm -f input.1 input.2 input.3 expected result 


echo "->goalign append phylip"
cat > input.1 <<EOF
   3   12
s1  ACGACGACGA CC
2  ATCTT-TTTT TC
3  ATCTT-TTTT TT
EOF

cat > input.2 <<EOF
   3   12
s4  ACGACGACGA CC
5  ATCTT-TTTT TC
6  ATCTT-TTTT TT
   3   12
s7  ACGACGACGA CC
8  ATCTT-TTTT TC
9  ATCTT-TTTT TT
EOF

cat > expected <<EOF
   9   12
s1  ACGACGACGA CC
2  ATCTT-TTTT TC
3  ATCTT-TTTT TT
s4  ACGACGACGA CC
5  ATCTT-TTTT TC
6  ATCTT-TTTT TT
s7  ACGACGACGA CC
8  ATCTT-TTTT TC
9  ATCTT-TTTT TT
EOF

${GOALIGN} append -i input.1 input.2 -p -o result
diff -q -b expected result
rm -f input.1 input.2 input.3 expected result 


echo "->goalign append phylip /2"
cat > input.1 <<EOF
   3   12
s1  ACGACGACGA CC
2  ATCTT-TTTT TC
3  ATCTT-TTTT TT
   3   12
s4  ACGACGACGA CC
5  ATCTT-TTTT TC
6  ATCTT-TTTT TT
   3   12
s7  ACGACGACGA CC
8  ATCTT-TTTT TC
9  ATCTT-TTTT TT
EOF

cat > expected <<EOF
   9   12
s1  ACGACGACGA CC
2  ATCTT-TTTT TC
3  ATCTT-TTTT TT
s4  ACGACGACGA CC
5  ATCTT-TTTT TC
6  ATCTT-TTTT TT
s7  ACGACGACGA CC
8  ATCTT-TTTT TC
9  ATCTT-TTTT TT
EOF

${GOALIGN} append -i input.1 -p -o result
diff -q -b expected result
rm -f input.1 input.2 input.3 expected result 

echo "->goalign stats char --per-sites"
cat > input <<EOF
>s1
ACGTACGT
>s2
ACGTACGT
>s3
ACGTACGT
>s4
ACGTTCGA
>s5
ACGTTCGA
EOF

cat > expected <<EOF
site	A	C	G	T
0	5	0	0	0
1	0	5	0	0
2	0	0	5	0
3	0	0	0	5
4	3	0	0	2
5	0	5	0	0
6	0	0	5	0
7	2	0	0	3
EOF

${GOALIGN} stats char --per-sites -i input > result
diff -q -b expected result
rm -f input expected result 

echo "->goalign stats char --per-sequences"
cat > input <<EOF
>s1
ACGTACGT
>s2
ACGTACGT
>s3
ACGTACGT
>s4
ACGTTCGA
>s5
ACGTTCGA
EOF

cat > expected <<EOF
seq	A	C	G	T
s1	2	2	2	2
s2	2	2	2	2
s3	2	2	2	2
s4	2	2	2	2
s5	2	2	2	2
EOF

${GOALIGN} stats char --per-sequences -i input > result
diff -q -b expected result
rm -f input expected result 



echo "->goalign stats char --per-sequences --count-profile"
cat > input <<EOF
>s1
ACGTACGT
>s2
-CGTACGT
>s3
ACGTACGT
>s4
ACGTTCGA
>s5
ACGTTCGA
EOF

cat > prof <<EOF
site	-	T	C	G	A
0	1	10	9	0	0
1	0	0	0	10	10
2	0	10	10	0	0
3	10	0	0	0	10
4	0	10	0	10	0
5	10	10	0	0	0
6	0	0	0	10	10
7	10	0	0	0	10
EOF

cat > expected <<EOF
sequence	gaps	gapsstart	gapsend	gapsuniques	gapsnew	gapsboth	gapsopenning	mutuniques	mutsnew	mutsboth	length	-	A	C	G	T
s1	0	0	0	0	0	0	0	0	7	0	8	0	2	2	2	2
s2	1	1	0	1	0	0	1	0	6	0	7	1	1	2	2	2
s3	0	0	0	0	0	0	0	0	7	0	8	0	2	2	2	2
s4	0	0	0	0	0	0	0	0	5	0	8	0	2	2	2	2
s5	0	0	0	0	0	0	0	0	5	0	8	0	2	2	2	2
EOF

${GOALIGN} stats --per-sequences -i input --count-profile prof  > result
diff -q -b expected result
rm -f input expected result 


echo "->goalign subsites /1"
cat > input <<EOF
>s1
ACGTACGT
>s2
-CGT-C-T
>s3
ACGTACGT
>s4
ACGTTCGA
>s5
ACGTTCGA
EOF

cat > exp <<EOF
>s1
CACGT
>s2
C-C-T
>s3
CACGT
>s4
CTCGA
>s5
CTCGA
EOF


${GOALIGN} subsites -i input 1 4 5 6 7 -o result
diff -q -b exp result
rm -f input exp result 

echo "->goalign subsites --reverse /2"
cat > input <<EOF
>s1
ACGTACGT
>s2
-CGT-C-T
>s3
ACGTACGT
>s4
ACGTTCGA
>s5
ACGTTCGA
EOF

cat > exp <<EOF
>s1
AGT
>s2
-GT
>s3
AGT
>s4
AGT
>s5
AGT
EOF

${GOALIGN} subsites -i input --reverse -o result 1 4 5 6 7 
diff -q -b exp result
rm -f input exp result 

echo "->goalign subsites --reverse --ref-seq /2"
cat > input <<EOF
>s1
ACGTACGT
>s2
-CGT-C-T
>s3
ACGTACGT
>s4
ACGTTCGA
>s5
ACGTTCGA
EOF

cat > exp <<EOF
>s1
ACTACG
>s2
-CT-C-
>s3
ACTACG
>s4
ACTTCG
>s5
ACTTCG
EOF

${GOALIGN} subsites -i input --reverse --ref-seq s2 -o result 1 4
diff -q -b exp result
rm -f input exp result 

echo "->goalign subsites from file"
cat > input <<EOF
>s1
ACGTACGT
>s2
-CGT-C-T
>s3
ACGTACGT
>s4
ACGTTCGA
>s5
ACGTTCGA
EOF

cat > posfile <<EOF
1
4
5
6
7
EOF


cat > exp <<EOF
>s1
CACGT
>s2
C-C-T
>s3
CACGT
>s4
CTCGA
>s5
CTCGA
EOF

${GOALIGN} subsites -i input --sitefile posfile -o result
diff -q -b exp result
rm -f input exp result 
