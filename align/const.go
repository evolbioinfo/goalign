package align

const (
	AMINOACIDS = 0 // Amino acid sequence alphabet
	NUCLEOTIDS = 1 // Nucleotid sequence alphabet
	UNKNOWN    = 2 // Unkown alphabet

	GAP   = '-'
	POINT = '.'
	OTHER = '*'

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
