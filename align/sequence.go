package align

import (
	"bytes"
	"errors"
	"fmt"
	"math/bits"
	"math/rand"
	"regexp"
	"strings"
	"unicode"
)

type Sequence interface {
	Sequence() string
	SequenceChar() []uint8
	SameSequence([]uint8) bool
	CharAt(int) uint8
	Name() string
	SetName(name string)
	Comment() string
	Length() int
	LongestORF() (start, end int) // Detects the longest ORF in forward strand only
	Reverse()
	Complement() error                                      // Returns an error if not nucleotide sequence
	Translate(phase int, geneticcode int) (Sequence, error) // Translates the sequence using the given code
	DetectAlphabet() int                                    // Try to detect alphabet (nt or aa)
	NumGaps() int                                           // Number of Gaps
	NumGapsOpenning() int                                   // Number of Gaps opennin, it counts streches of gap only once
	NumGapsFromStart() int                                  // Number of Gaps from Start (until a non gap is encountered)
	NumGapsFromEnd() int                                    // Number of Gaps from End (until a non gap is encountered)
	// returns the number of differences between the reference sequence and each sequence of the alignment
	// If lengths are different, returns an error
	// It does not take into account 'N' and '-' in sequences as mutations compared to ref
	/// sequence (ref sequence can have a '-' or a 'N')
	NumMutationsComparedToReferenceSequence(alphabet int, seq Sequence) (nummutations int, err error)
	// returns the list ofdifferences between the reference sequence and each sequence of the alignment
	// Counts only non N sites in each sequences (may be a gap or a N in the reference sequence though)
	// If a character is ambigous (IUPAC notation), then it is counted as a mutation only if it is incompatible with
	// the reference character.
	// if aa is true: the sequences are nucleotides and nucleotides are taken codon by codon of the reference sequence
	// to list mutations. In case of insertion or a deletion in the target sequence: if %3==0: - or aa insert,
	// otherwise "/" ~frameshift?
	//
	// If lengths are different, returns an error
	ListMutationsComparedToReferenceSequence(alphabet int, refseq Sequence, aa bool) (mutations []Mutation, err error)

	Clone() Sequence
}

type seq struct {
	name     string  // Name of the sequence
	sequence []uint8 // Sequence of nucleotides/aa
	comment  string  // Comment if any
}

func NewSequence(name string, sequence []uint8, comment string) *seq {
	return &seq{
		name,
		sequence,
		comment,
	}
}

func (s *seq) Sequence() string {
	return string(s.sequence)
}
func (s *seq) SequenceChar() []uint8 {
	return s.sequence
}

/*
Returns true iif :
- slices have the same length and
- slices have same uint8 at each position
*/
func (s *seq) SameSequence(runeseq []uint8) bool {
	if len(s.sequence) != len(runeseq) {
		return false
	}
	for i, v := range s.sequence {
		if runeseq[i] != v {
			return false
		}
	}
	return true
}

func (s *seq) CharAt(i int) uint8 {
	return s.sequence[i]
}

func (s *seq) Name() string {
	return s.name
}

func (s *seq) SetName(name string) {
	s.name = name
}

func (s *seq) Comment() string {
	return s.comment
}

func (s *seq) Length() int {
	return len(s.sequence)
}

// Detects the position of ATG giving the longest ORF
// Search is done in the forward strand only
//
// returns -1 is no ATG...STOP has been found
func (s *seq) LongestORF() (start, end int) {
	start = -1
	end = -1
	re, _ := regexp.Compile("(ATG)(.{3})*?(TAA|TGA|TAG)")
	//re.Longest()
	idx := re.FindAllStringIndex(
		strings.Replace(
			strings.ToUpper(string(s.sequence)),
			"U", "T", -1),
		-1)
	for _, pos := range idx {
		if pos[1]-pos[0] > end-start {
			end = pos[1]
			start = pos[0]
		}
	}
	return start, end
}

func RandomSequence(alphabet, length int) ([]uint8, error) {
	seq := make([]uint8, length)
	for i := 0; i < length; i++ {
		switch alphabet {
		case AMINOACIDS:
			seq[i] = stdaminoacid[rand.Intn(len(stdaminoacid))]
		case NUCLEOTIDS:
			seq[i] = stdnucleotides[rand.Intn(len(stdnucleotides))]
		default:
			return nil, errors.New("unexpected sequence alphabet type")
		}
	}
	return seq, nil
}

// Reverses a sequence
func Reverse(seq []uint8) {
	for i, j := 0, len(seq)-1; i < j; i, j = i+1, j-1 {
		seq[i], seq[j] = seq[j], seq[i]
	}
}

// Reverse sequence order
func (s *seq) Reverse() {
	Reverse(s.sequence)
}

// Complement sequence
func Complement(seq []uint8) (err error) {

	for i, n := range seq {
		c, ok := complement_nuc_mapping[n]
		if !ok {
			err = fmt.Errorf("character %c can not be complemented", n)
			return
		}
		seq[i] = c
	}
	return
}

// Complement sequence
func (s *seq) Complement() error {
	a := s.DetectAlphabet()
	if a != NUCLEOTIDS && a != BOTH {
		return fmt.Errorf("wrong alphabet for Complementing sequence")
	}
	return Complement(s.sequence)
}

func (s *seq) DetectAlphabet() int {
	isaa := true
	isnt := true

	for _, nt := range s.sequence {
		nt = uint8(unicode.ToUpper(rune(nt)))
		couldbent := false
		couldbeaa := false
		switch nt {
		case 'A', 'C', 'B', 'R', 'G', '?', GAP, POINT, OTHER, 'D', 'K', 'S', 'H', 'M', 'N', 'V', 'X', 'T', 'W', 'Y':
			couldbent = true
			couldbeaa = true
		case 'U', 'O':
			couldbent = true
		case 'Q', 'E', 'I', 'L', 'F', 'P', 'Z':
			couldbeaa = true
		}
		isaa = isaa && couldbeaa
		isnt = isnt && couldbent
	}

	if isnt && isaa {
		return BOTH
	} else if isnt {
		return NUCLEOTIDS
	} else if isaa {
		return AMINOACIDS
	} else {
		return UNKNOWN
	}
}

// NumGaps returns the number of Gaps on the given sequence
func (s *seq) NumGaps() (numgaps int) {
	numgaps = 0
	for _, c := range s.sequence {
		if c == GAP {
			numgaps++
		}
	}
	return
}

// NumGapsOpenning returns the number of Gaps on the given sequence
func (s *seq) NumGapsOpenning() (numgaps int) {
	numgaps = 0
	var prevChar uint8 = '>'
	for _, c := range s.sequence {
		if c == GAP && prevChar != GAP {
			numgaps++
		}
		prevChar = c
	}
	return
}

// NumGapsFromStart returns the number of Gaps on the given sequence
// by counting only the first consecutive gaps
// Ex: -----A-AAA--AA = 5
func (s *seq) NumGapsFromStart() (numgaps int) {
	numgaps = 0
	for _, c := range s.sequence {
		if c != GAP {
			return
		}
		numgaps++
	}
	return
}

// NumGapsFromEnd returns the number of Gaps on the given sequence
// by counting only the last consecutive gaps at the end
// Ex: // -----A-AAA--AA---- = 4
func (s *seq) NumGapsFromEnd() (numgaps int) {
	numgaps = 0
	for i := range s.sequence {
		c := s.sequence[len(s.sequence)-i-1]
		if c != GAP {
			return
		}
		numgaps++
	}
	return
}

// returns the number of differences between the reference sequence and each sequence of the alignment
// Counts only non GAPS and non N sites in each sequences (may be a gap or a N in the reference sequence though)
// If a character is ambigous (IUPAC notation), then it is counted as a mutation only if it is incompatible with
// the reference character.
//
// If lengths are different, returns an error
func (s *seq) NumMutationsComparedToReferenceSequence(alphabet int, refseq Sequence) (nummutations int, err error) {
	var refseqCode []uint8
	var nt uint8

	if refseq.Length() != s.Length() {
		err = fmt.Errorf("reference sequence and sequence do not have same length (%d,%d), cannot compute a number of mutation", refseq.Length(), s.Length())
		return
	}

	all := uint8('.')
	if alphabet == NUCLEOTIDS {
		refseqCode = make([]uint8, s.Length())
		for i := 0; i < s.Length(); i++ {
			if refseqCode[i], err = Nt2IndexIUPAC(refseq.SequenceChar()[i]); err != nil {
				return
			}
		}
		all = ALL_NUCLE
	} else {
		all = ALL_AMINO
	}

	for i := 0; i < s.Length(); i++ {
		eq := true
		if alphabet == NUCLEOTIDS {
			if nt, err = Nt2IndexIUPAC(s.sequence[i]); err != nil {
				return
			}
			if eq, err = EqualOrCompatible(nt, refseqCode[i]); err != nil {
				return
			}
		} else {
			eq = (s.sequence[i] == refseq.SequenceChar()[i])
		}
		if s.SequenceChar()[i] != GAP && s.SequenceChar()[i] != all && !eq {
			nummutations++
		}
	}
	return
}

// returns the list ofdifferences between the reference sequence and each sequence of the alignment
// Counts only non N sites in each sequences (may be a gap or a N in the reference sequence though)
// If a character is ambigous (IUPAC notation), then it is counted as a mutation only if it is incompatible with
// the reference character.
//
// If lengths are different, returns an error
func (s *seq) ListMutationsComparedToReferenceSequence(alphabet int, refseq Sequence, aa bool) (mutations []Mutation, err error) {
	if aa {
		return s.listMutationsComparedToReferenceSequenceAA(alphabet, refseq)
	}
	return s.listMutationsComparedToReferenceSequence(alphabet, refseq)
}

func (s *seq) listMutationsComparedToReferenceSequence(alphabet int, refseq Sequence) (mutations []Mutation, err error) {
	var refseqCode []uint8
	var refseqchar []uint8
	var nt uint8

	if refseq.Length() != s.Length() {
		err = fmt.Errorf("reference sequence and sequence do not have same length (%d,%d), cannot compute a number of mutation", refseq.Length(), s.Length())
		return
	}

	refseqchar = refseq.SequenceChar()
	all := uint8('.')
	if alphabet == NUCLEOTIDS {
		refseqCode = make([]uint8, s.Length())
		for i := 0; i < s.Length(); i++ {
			if refseqCode[i], err = Nt2IndexIUPAC(refseqchar[i]); err != nil {
				return
			}
		}
		all = ALL_NUCLE
	} else {
		all = ALL_AMINO
	}

	curinsert := make([]uint8, 0) // buffer of current consecutive insertions in query seq (gaps in ref seq)
	refi := 0                     // Current index on ref sequence (without gaps)
	for i := 0; i < s.Length(); i++ {
		eq := true
		if alphabet == NUCLEOTIDS {
			if nt, err = Nt2IndexIUPAC(s.sequence[i]); err != nil {
				return
			}
			if eq, err = EqualOrCompatible(nt, refseqCode[i]); err != nil {
				return
			}
		} else {
			eq = (s.sequence[i] == refseqchar[i])
		}

		if refseqchar[i] == GAP {
			if s.sequence[i] != GAP {
				curinsert = append(curinsert, s.sequence[i])
			}
		} else {
			if len(curinsert) > 0 {
				mutations = append(mutations, Mutation{Ref: '-', Pos: refi, Alt: curinsert})
				curinsert = make([]uint8, 0)
			}
			if s.sequence[i] != all && !eq {
				mutations = append(mutations, Mutation{Ref: refseqchar[i], Pos: refi, Alt: []uint8{s.sequence[i]}})
			}
		}

		if refseqchar[i] != GAP {
			refi++
		}
	}
	if len(curinsert) > 0 {
		mutations = append(mutations, Mutation{Ref: '-', Pos: refi, Alt: curinsert})
		curinsert = make([]uint8, 0)
	}

	return
}

// listMutationsComparedToReferenceSequenceAA takes nucleotides codon by codon of the reference seqence
// to list mutations. In case of insertion or a deletion in the target sequence: if %3==0: - or aa insert,
// otherwise "/" ~frameshift?
func (s *seq) listMutationsComparedToReferenceSequenceAA(alphabet int, refseq Sequence) (mutations []Mutation, err error) {
	var refseqchar []uint8
	var code map[string]uint8

	if refseq.Length() != s.Length() {
		err = fmt.Errorf("reference sequence and sequence do not have same length (%d,%d), cannot compute a number of mutation", refseq.Length(), s.Length())
		return
	}

	if alphabet != NUCLEOTIDS {
		err = fmt.Errorf("alignment must be nucleotidic to work codon by codon")
		return
	}

	if code, err = geneticCode(GENETIC_CODE_STANDARD); err != nil {
		return
	}

	refseqchar = refseq.SequenceChar()

	refcodonidx := []int{0, 1, 2}
	// We traverse reference nt 3 by 3
	//The reference codon may have gaps between nt ,
	// ex 1:
	// Ref: AC--GTACGT
	// Seq: ACTTGTACGT
	// In that case, the first ref codon is [0,1,4], corresponding to sequence ACTTG in seq
	// ACTTG % 3 != 0 ==> Frameshift?
	// ex 2:
	// Ref: AC---GTACGT
	// Seq: ACTTTGTACGT
	// ref codon: [0,1,5]
	// seq      : ACTTTG : Insertion - OK
	// ex 3:
	// Ref: ACGTACGT
	// Seq: A--TACGT
	// ref codon: [0,1,2]
	// ses      : A--: Deletion: not ok : Frameshift?
	aaidx := 0
	for refcodonidx[2] < s.Length() {
		var refaa uint8
		var allgaps bool = false

		if refseqchar[refcodonidx[0]] == GAP && refseqchar[refcodonidx[1]] == GAP && refseqchar[refcodonidx[2]] == GAP {
			refaa = '-'
			aaidx--
			allgaps = true
		} else {
			// We find the three reference positions without gap
			for refcodonidx[2] < s.Length() && refseqchar[refcodonidx[0]] == GAP {
				refcodonidx[0]++
				refcodonidx[1]++
				refcodonidx[2]++
			}
			if refcodonidx[2] >= s.Length() {
				break
			}
			for refcodonidx[2] < s.Length() && refseqchar[refcodonidx[1]] == GAP {
				refcodonidx[1]++
				refcodonidx[2]++
			}
			if refcodonidx[2] >= s.Length() {
				break
			}
			for refcodonidx[2] < s.Length() && refseqchar[refcodonidx[2]] == GAP {
				refcodonidx[2]++
			}
			if refcodonidx[2] >= s.Length() {
				break
			}
			refaa = translateCodon(refseqchar[refcodonidx[0]], refseqchar[refcodonidx[1]], refseqchar[refcodonidx[2]], code)
		}
		// We remove the gaps from the compared sequence corresponding to the reference codon
		tmpseq := make([]uint8, 0)
		for si := refcodonidx[0]; si <= refcodonidx[2]; si++ {
			if s.sequence[si] != GAP {
				tmpseq = append(tmpseq, s.sequence[si])
			}
		}

		// Deletion
		if len(tmpseq) == 0 {
			if !allgaps {
				mutations = append(mutations, Mutation{Ref: refaa, Pos: aaidx, Alt: []uint8{'-'}})
			}
		} else if len(tmpseq)%3 != 0 {
			// Potential frameshift
			mutations = append(mutations, Mutation{Ref: refaa, Pos: aaidx, Alt: []uint8{'/'}})
		} else {
			// Mutations + potentiel Insertions
			// We find all corresponding codons of the target sequence, between refcodonidx[0] and refcodonidx[2]
			curinsert := make([]uint8, 0) // buffer of current consecutive insertions in query seq (gaps in ref seq)
			diff := false
			for si := 0; si <= len(tmpseq)-2; si += 3 {
				aa := translateCodon(tmpseq[si], tmpseq[si+1], tmpseq[si+2], code)
				curinsert = append(curinsert, aa)
				diff = diff || aa != refaa
			}
			if len(curinsert) > 1 || diff {
				mutations = append(mutations, Mutation{Ref: refaa, Pos: aaidx, Alt: curinsert})
			}
		}

		refcodonidx[0] = refcodonidx[2] + 1
		refcodonidx[1] = refcodonidx[2] + 2
		refcodonidx[2] = refcodonidx[2] + 3
		aaidx++
	}
	return
}

// Translates the given sequence start at nucleotide index "phase"
//
// If the sequence is not nucleotidic, then throws an error
// If sequence length is < 3+phase then thropws an error
func (s *seq) Translate(phase int, geneticcode int) (tr Sequence, err error) {
	var buffer bytes.Buffer
	var code map[string]uint8

	if code, err = geneticCode(geneticcode); err != nil {
		return
	}

	if err = bufferTranslate(s, phase, code, &buffer); err != nil {
		return
	}

	tr = NewSequence(s.name, []uint8(buffer.String()), s.comment)
	return
}

func bufferTranslate(s *seq, phase int, code map[string]uint8, buffer *bytes.Buffer) (err error) {
	buffer.Reset()
	if s.DetectAlphabet() != NUCLEOTIDS && s.DetectAlphabet() != BOTH {
		err = fmt.Errorf("cannot translate this sequence, wrong alphabet")
		return
	}

	if len(s.sequence) < 3+phase {
		err = fmt.Errorf("cannot translate a sequence with length < 3+phase (%s)", s.name)
		return
	}
	for i := phase; i < len(s.sequence)-2; i += 3 {
		var aa uint8 = ' '
		aa = translateCodon(s.sequence[i], s.sequence[i+1], s.sequence[i+2], code)
		buffer.WriteByte(aa)
	}
	return
}

// We translate the given codon
// If IUPAC codes : we generate all possible codons.
// If all codons give the same amino acid: we take this one.
// Otherwise we translate by 'X'
func translateCodon(n1, n2, n3 uint8, code map[string]uint8) (aa uint8) {
	var aatmp uint8 = ' '
	var found bool = false
	aa = ' '
	// We handle possible IUPAC characters
	codons := GenAllPossibleCodons(n1, n2, n3)
	if len(codons) == 0 {
		aa = 'X'
	}
	for _, codon := range codons {
		// The codon i=s not found
		// We return X
		if aatmp, found = code[codon]; !found {
			aa = 'X'
			break
		}
		// Different codons give different AA
		// We can not translate it uniquely
		if aa != ' ' && aatmp != aa {
			aa = 'X'
			break
		}
		aa = aatmp
	}
	return
}

func (s *seq) Clone() Sequence {
	seq2 := make([]uint8, len(s.sequence))
	copy(seq2, s.sequence)
	return NewSequence(s.name, seq2, s.comment)
}

// GenAllPossibleCodons generates all possible codons given the 3 nucleotides in arguments
// Multiple codons may exist if IUPAC code is employed (R=A|G, etc.).
// The 3 nucleotites in arguments are converted to upper case and U converted to T.
// If one character does not correspond to a known nucleotide in IUPAC code, then
// Returns an empty slice.
// If one of the nucleotides is a GAP, then returns an empty slice.
//
// For example GenAllPossibleCodons('A','G','N') should return
// {"AGA","AGC","AGG","AGT"}.
func GenAllPossibleCodons(nt1, nt2, nt3 uint8) (codons []string) {
	var nt uint8
	var codontmp string
	var nts1, nts2, nts3 []uint8 // possible nts for each nt
	var found bool

	codons = make([]string, 0)
	codonstmp := make([]string, 0)

	nt1 = uint8(unicode.ToUpper(rune(nt1)))
	nt2 = uint8(unicode.ToUpper(rune(nt2)))
	nt3 = uint8(unicode.ToUpper(rune(nt3)))

	if nt1 == 'U' {
		nt1 = 'T'
	}
	if nt2 == 'U' {
		nt2 = 'T'
	}
	if nt3 == 'U' {
		nt3 = 'T'
	}

	if nts1, found = IupacCode[nt1]; !found {
		return
	}
	if nts2, found = IupacCode[nt2]; !found {
		return
	}
	if nts3, found = IupacCode[nt3]; !found {
		return
	}

	for _, nt = range nts1 {
		codons = append(codons, string(nt))
	}
	for _, nt = range nts2 {
		for _, codontmp = range codons {
			codontmp = codontmp + string(nt)
			codonstmp = append(codonstmp, string(codontmp))
		}
	}
	codons = codonstmp
	codonstmp = make([]string, 0)

	for _, nt = range nts3 {
		for _, codontmp = range codons {
			codontmp = codontmp + string(nt)
			codonstmp = append(codonstmp, string(codontmp))
		}
	}
	codons = codonstmp
	return
}

// EqualOrCompatible returns true if the two nucleotides are identical or
// if they are compatible in case they are ambigous.
//
// For example :
// Y: {C | T} is compatible with S: {G | C} because there is one nt in common
// If nt1 or nt2 are not nucleotides, then returns an error
// n1 and nt2 valures are from NT_... in const.go
func EqualOrCompatible(nt1, nt2 uint8) (ok bool, err error) {

	if nt1 > NT_N {
		err = fmt.Errorf("given nucleotide 1 code (%d) does not exist", nt1)
		return
	}
	if nt2 > NT_N {
		err = fmt.Errorf("given nucleotide 1 code (%d) does not exist", nt1)
		return
	}
	if nt1 == nt2 {
		ok = true
		return
	}
	ok = (nt1 & nt2) > 0

	return
}

// NtIUPACDifference returns the cost of the difference between
// the two potentially ambiguous nucleotides.
//
// - if the two nucleotides are identical : returns 0.0
// - if the two nucleotides are different:
//  1. If none are ambigous: returns 1.0
//  2. Otherwise, returns 1-Card(I)/Card(U), I being the
//     intersection of the sets of possible
//     nucleotides of nt1 and nt2, and U being
//     the union of the sets of possible nucleotides
//     of nt1 and nt2.
//
// For example, if we want to compare Y and S :
// Y = {C | T} and S = {G | C}. Card(I)=1, Card(U)=3, so diff=2/3
//
// Precisions:
// - For N vs. A for example: the difference will be 1-1/4 : 3/4
// - For gaps: Returns diff=1.0
//
// nt1 and nt2 values are in NT_... of const.go
func NtIUPACDifference(nt1, nt2 uint8) (diff float64, err error) {
	diff = 0.0

	if nt1 > NT_N {
		err = fmt.Errorf("given nucleotide 1 code (%d) does not exist", nt1)
		return
	}
	if nt2 > NT_N {
		err = fmt.Errorf("given nucleotide 1 code (%d) does not exist", nt1)
		return
	}

	if nt1 == nt2 {
		return
	}
	inter := bits.OnesCount8(nt1 & nt2)
	//union := bits.OnesCount8(nt1 | nt2)
	// diff = 1.0 - float64(inter)/float64(union)
	if inter == 0 {
		diff = 1.0
	}
	return
}
