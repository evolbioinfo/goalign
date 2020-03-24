package align

import (
	"errors"
	"fmt"
	"unicode"
)

const (
	AMINOACIDS = 0 // Amino acid sequence alphabet
	NUCLEOTIDS = 1 // Nucleotid sequence alphabet
	BOTH       = 2 // Could be both
	UNKNOWN    = 3 // Unkown alphabet

	GAP   = '-'
	POINT = '.'
	OTHER = '*'

	ALL_AMINO = 'X'
	ALL_NUCLE = 'N'

	PSSM_NORM_NONE = 0 // No normalization
	PSSM_NORM_FREQ = 1 // Normalization by freq in the site
	PSSM_NORM_DATA = 2 // Normalization by aa/nt frequency in data
	PSSM_NORM_UNIF = 3 // Normalization by uniform frequency
	PSSM_NORM_LOGO = 4 // Normalization like LOGO : v(site)=freq*(log2(alphabet)-H(site)-pseudocount

	FORMAT_FASTA   = 0
	FORMAT_PHYLIP  = 1
	FORMAT_NEXUS   = 2
	FORMAT_CLUSTAL = 3

	POSITION_IDENTICAL      = 0 // All characters in a position are the same
	POSITION_CONSERVED      = 1 // Same strong group
	POSITION_SEMI_CONSERVED = 2 // Same weak group
	POSITION_NOT_CONSERVED  = 3 // None of the above values

	GENETIC_CODE_STANDARD         = 0 // Standard genetic code
	GENETIC_CODE_VETEBRATE_MITO   = 1 // Vertebrate mitochondrial genetic code
	GENETIC_CODE_INVETEBRATE_MITO = 2 // Invertebrate mitochondrial genetic code
)

var stdaminoacid = []rune{'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'}
var stdnucleotides = []rune{'A', 'C', 'G', 'T'}

var standardcode = map[string]rune{
	"---": '-',
	"GCT": 'A', "GCC": 'A', "GCA": 'A', "GCG": 'A',
	"TTA": 'L', "TTG": 'L', "CTT": 'L', "CTC": 'L', "CTA": 'L', "CTG": 'L',
	"CGT": 'R', "CGC": 'R', "CGA": 'R', "CGG": 'R', "AGA": 'R', "AGG": 'R',
	"AAA": 'K', "AAG": 'K',
	"AAT": 'N', "AAC": 'N',
	"ATG": 'M',
	"GAT": 'D', "GAC": 'D',
	"TTT": 'F', "TTC": 'F',
	"TGT": 'C', "TGC": 'C',
	"CCT": 'P', "CCC": 'P', "CCA": 'P', "CCG": 'P',
	"CAA": 'Q', "CAG": 'Q',
	"TCT": 'S', "TCC": 'S', "TCA": 'S', "TCG": 'S', "AGT": 'S', "AGC": 'S',
	"GAA": 'E', "GAG": 'E',
	"ACT": 'T', "ACC": 'T', "ACA": 'T', "ACG": 'T',
	"GGT": 'G', "GGC": 'G', "GGA": 'G', "GGG": 'G',
	"TGG": 'W',
	"CAT": 'H', "CAC": 'H',
	"TAT": 'Y', "TAC": 'Y',
	"ATT": 'I', "ATC": 'I', "ATA": 'I',
	"GTT": 'V', "GTC": 'V', "GTA": 'V', "GTG": 'V',
	"TAA": '*', "TGA": '*', "TAG": '*',
}

var vertebratemitocode = map[string]rune{
	"---": '-',
	"GCT": 'A', "GCC": 'A', "GCA": 'A', "GCG": 'A',
	"TTA": 'L', "TTG": 'L', "CTT": 'L', "CTC": 'L', "CTA": 'L', "CTG": 'L',
	"CGT": 'R', "CGC": 'R', "CGA": 'R', "CGG": 'R',
	"AAA": 'K', "AAG": 'K',
	"AAT": 'N', "AAC": 'N',
	"ATG": 'M', "ATA": 'M',
	"GAT": 'D', "GAC": 'D',
	"TTT": 'F', "TTC": 'F',
	"TGT": 'C', "TGC": 'C',
	"CCT": 'P', "CCC": 'P', "CCA": 'P', "CCG": 'P',
	"CAA": 'Q', "CAG": 'Q',
	"TCT": 'S', "TCC": 'S', "TCA": 'S', "TCG": 'S', "AGT": 'S', "AGC": 'S',
	"GAA": 'E', "GAG": 'E',
	"ACT": 'T', "ACC": 'T', "ACA": 'T', "ACG": 'T',
	"GGT": 'G', "GGC": 'G', "GGA": 'G', "GGG": 'G',
	"TGG": 'W', "TGA": 'W',
	"CAT": 'H', "CAC": 'H',
	"TAT": 'Y', "TAC": 'Y',
	"ATT": 'I', "ATC": 'I',
	"GTT": 'V', "GTC": 'V', "GTA": 'V', "GTG": 'V',
	"AGA": '*', "AGG": '*', "TAA": '*', "TAG": '*',
}

var invertebratemitocode = map[string]rune{
	"---": '-',
	"GCT": 'A', "GCC": 'A', "GCA": 'A', "GCG": 'A',
	"TTA": 'L', "TTG": 'L', "CTT": 'L', "CTC": 'L', "CTA": 'L', "CTG": 'L',
	"CGT": 'R', "CGC": 'R', "CGA": 'R', "CGG": 'R',
	"AAA": 'K', "AAG": 'K',
	"AAT": 'N', "AAC": 'N',
	"ATG": 'M', "ATA": 'M',
	"GAT": 'D', "GAC": 'D',
	"TTT": 'F', "TTC": 'F',
	"TGT": 'C', "TGC": 'C',
	"CCT": 'P', "CCC": 'P', "CCA": 'P', "CCG": 'P',
	"CAA": 'Q', "CAG": 'Q',
	"AGA": 'S', "AGG": 'S', "TCT": 'S', "TCC": 'S', "TCA": 'S', "TCG": 'S', "AGT": 'S', "AGC": 'S',
	"GAA": 'E', "GAG": 'E',
	"ACT": 'T', "ACC": 'T', "ACA": 'T', "ACG": 'T',
	"GGT": 'G', "GGC": 'G', "GGA": 'G', "GGG": 'G',
	"TGG": 'W', "TGA": 'W',
	"CAT": 'H', "CAC": 'H',
	"TAT": 'Y', "TAC": 'Y',
	"ATT": 'I', "ATC": 'I',
	"GTT": 'V', "GTC": 'V', "GTA": 'V', "GTG": 'V',
	"TAA": '*', "TAG": '*',
}

func geneticCode(code int) (gencode map[string]rune, err error) {
	switch code {
	case GENETIC_CODE_STANDARD:
		gencode = standardcode
	case GENETIC_CODE_VETEBRATE_MITO:
		gencode = vertebratemitocode
	case GENETIC_CODE_INVETEBRATE_MITO:
		gencode = invertebratemitocode
	default:
		err = fmt.Errorf("This genetic code does not exis")
	}
	return
}

// Amino Acid strong groups for clustal format
// conservation line
//
// Taken from Clustalw source code
var strongGroups = [][]rune{
	{'S', 'T', 'A'},
	{'N', 'E', 'Q', 'K'},
	{'N', 'H', 'Q', 'K'},
	{'N', 'D', 'E', 'Q'},
	{'Q', 'H', 'R', 'K'},
	{'M', 'I', 'L', 'V'},
	{'M', 'I', 'L', 'F'},
	{'H', 'Y'},
	{'F', 'Y', 'W'},
}

// Amino Acid weak groups for clustal format
// conservation line
//
// Taken from Clustalw source code
var weakGroups = [][]rune{
	{'C', 'S', 'A'},
	{'A', 'T', 'V'},
	{'S', 'A', 'G'},
	{'S', 'T', 'N', 'K'},
	{'S', 'T', 'P', 'A'},
	{'S', 'G', 'N', 'D'},
	{'S', 'N', 'D', 'E', 'Q', 'K'},
	{'N', 'D', 'E', 'Q', 'H', 'K'},
	{'N', 'E', 'Q', 'H', 'R', 'K'},
	{'F', 'V', 'L', 'I', 'M'},
	{'H', 'F', 'Y'},
}

var dna_to_matrix_pos = map[rune]int{
	'A': 0,
	'T': 1,
	'G': 2,
	'C': 3,
	'S': 4,
	'W': 5,
	'R': 6,
	'Y': 7,
	'K': 8,
	'M': 9,
	'B': 10,
	'V': 11,
	'H': 12,
	'D': 13,
	'N': 14,
	'U': 15,
	'X': 14,
}

var iupacCode = map[rune][]rune{
	'A': {'A'},
	'C': {'C'},
	'G': {'G'},
	'T': {'T'},
	'R': {'A', 'G'},
	'Y': {'C', 'T'},
	'S': {'G', 'C'},
	'W': {'A', 'T'},
	'K': {'G', 'T'},
	'M': {'A', 'C'},
	'B': {'C', 'G', 'T'},
	'D': {'A', 'G', 'T'},
	'H': {'A', 'C', 'T'},
	'V': {'A', 'C', 'G'},
	'N': {'A', 'C', 'G', 'T'},
}

// Taken from EMBOSS Water
// This matrix was created by Todd Lowe   12/10/92
//
// Uses ambiguous nucleotide codes, probabilities rounded to
//  nearest integer
//
// Lowest score = -4, Highest score = 5
//
//A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N   U
var dnafull_subst_matrix = [][]float64{
	[]float64{5, -4, -4, -4, -4, 1, 1, -4, -4, 1, -4, -1, -1, -1, -2, -4},
	[]float64{-4, 5, -4, -4, -4, 1, -4, 1, 1, -4, -1, -4, -1, -1, -2, 5},
	[]float64{-4, -4, 5, -4, 1, -4, 1, -4, 1, -4, -1, -1, -4, -1, -2, -4},
	[]float64{-4, -4, -4, 5, 1, -4, -4, 1, -4, 1, -1, -1, -1, -4, -2, -4},
	[]float64{-4, -4, 1, 1, -1, -4, -2, -2, -2, -2, -1, -1, -3, -3, -1, -4},
	[]float64{1, 1, -4, -4, -4, -1, -2, -2, -2, -2, -3, -3, -1, -1, -1, 1},
	[]float64{1, -4, 1, -4, -2, -2, -1, -4, -2, -2, -3, -1, -3, -1, -1, -4},
	[]float64{-4, 1, -4, 1, -2, -2, -4, -1, -2, -2, -1, -3, -1, -3, -1, 1},
	[]float64{-4, 1, 1, -4, -2, -2, -2, -2, -1, -4, -1, -3, -3, -1, -1, 1},
	[]float64{1, -4, -4, 1, -2, -2, -2, -2, -4, -1, -3, -1, -1, -3, -1, -4},
	[]float64{-4, -1, -1, -1, -1, -3, -3, -1, -1, -3, -1, -2, -2, -2, -1, -1},
	[]float64{-1, -4, -1, -1, -1, -3, -1, -3, -3, -1, -2, -1, -2, -2, -1, -4},
	[]float64{-1, -1, -4, -1, -3, -1, -3, -1, -3, -1, -2, -2, -1, -2, -1, -1},
	[]float64{-1, -1, -1, -4, -3, -1, -1, -3, -1, -3, -2, -2, -2, -1, -1, -1},
	[]float64{-2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2},
	[]float64{-4, 5, -4, -4, -4, 1, -4, 1, 1, -4, -1, -4, -1, -1, -2, 5},
}

var prot_to_matrix_pos = map[rune]int{
	'A': 0,
	'R': 1,
	'N': 2,
	'D': 3,
	'C': 4,
	'Q': 5,
	'E': 6,
	'G': 7,
	'H': 8,
	'I': 9,
	'L': 10,
	'K': 11,
	'M': 12,
	'F': 13,
	'P': 14,
	'S': 15,
	'T': 16,
	'W': 17,
	'Y': 18,
	'V': 19,
	'B': 20,
	'Z': 21,
	'X': 22,
	'*': 23,
}

// Taken from EMBOSS WATER
// Matrix made by matblas from blosum62.iij
// * column uses minimum score
// BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
// Blocks Database = /data/blocks_5.0/blocks.dat
// Cluster Percentage: >= 62
// Entropy =   0.6979, Expected =  -0.5209
// A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
var blosum62_subst_matrix = [][]float64{
	[]float64{4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0, -4},
	[]float64{-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1, -4},
	[]float64{-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1, -4},
	[]float64{-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4},
	[]float64{0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
	[]float64{-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1, -4},
	[]float64{-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4},
	[]float64{0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1, -4},
	[]float64{-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4},
	[]float64{-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4},
	[]float64{-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -4},
	[]float64{-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1, -4},
	[]float64{-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4},
	[]float64{-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4},
	[]float64{-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
	[]float64{1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0, -4},
	[]float64{0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0, -4},
	[]float64{-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2, -4},
	[]float64{-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -4},
	[]float64{0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -4},
	[]float64{-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1, -4},
	[]float64{-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4},
	[]float64{0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4},
	[]float64{-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1},
}

var complement_nuc_mapping = map[rune]rune{
	'A': 'T',
	'T': 'A',
	'U': 'A',
	'G': 'C',
	'C': 'G',
	'Y': 'R',
	'R': 'Y',
	'S': 'S',
	'W': 'W',
	'K': 'M',
	'M': 'K',
	'B': 'V',
	'D': 'H',
	'H': 'D',
	'V': 'B',
	'N': 'N',
}

/*
Returns the index of each nts
A=0
C=1
G=2
T=3
*/
func Nt2Index(nt rune) (idx int, err error) {
	switch unicode.ToUpper(nt) {
	case 'A':
		idx = 0
	case 'C':
		idx = 1
	case 'G':
		idx = 2
	case 'T':
		idx = 3
	default:
		err = errors.New(fmt.Sprintf("No index for character: %c", nt))
	}
	return
}

/*
Returns the index of each nts
0=A
1=C
2=G
3=T
*/
func Index2Nt(index int) (nt rune, err error) {
	if index < 0 || index > len(stdnucleotides) {
		err = errors.New(fmt.Sprintf("No nt for index %d", index))
	}
	nt = stdnucleotides[index]
	return
}

func AA2Index(aa rune) (idx int, err error) {
	switch unicode.ToUpper(aa) {
	case 'A':
		idx = 0
	case 'R':
		idx = 1
	case 'N':
		idx = 2
	case 'D':
		idx = 3
	case 'C':
		idx = 4
	case 'Q':
		idx = 5
	case 'E':
		idx = 6
	case 'G':
		idx = 7
	case 'H':
		idx = 8
	case 'I':
		idx = 9
	case 'L':
		idx = 10
	case 'K':
		idx = 11
	case 'M':
		idx = 12
	case 'F':
		idx = 13
	case 'P':
		idx = 14
	case 'S':
		idx = 15
	case 'T':
		idx = 16
	case 'W':
		idx = 17
	case 'Y':
		idx = 18
	case 'V':
		idx = 19
	default:
		err = errors.New(fmt.Sprintf("No index for aa: %c", aa))
	}
	return
}

func Index2AA(index int) (aa rune, err error) {
	if index < 0 || index > len(stdaminoacid) {
		err = errors.New(fmt.Sprintf("No aa for index %d", index))
	}
	aa = stdaminoacid[index]
	return
}
