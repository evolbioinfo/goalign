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
	SequenceChar() []rune
	SameSequence([]rune) bool
	CharAt(int) rune
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
	NumGapsFromStart() int                                  // Number of Gaps from Start (until a non gap is encountered)
	NumGapsFromEnd() int                                    // Number of Gaps from End (until a non gap is encountered)
	Clone() Sequence
}

type seq struct {
	name     string // Name of the sequence
	sequence []rune // Sequence of nucleotides/aa
	comment  string // Comment if any
}

func NewSequence(name string, sequence []rune, comment string) *seq {
	return &seq{
		name,
		sequence,
		comment,
	}
}

func (s *seq) Sequence() string {
	return string(s.sequence)
}
func (s *seq) SequenceChar() []rune {
	return s.sequence
}

/*
Returns true iif :
- slices have the same length and
- slices have same rune at each position
*/
func (s *seq) SameSequence(runeseq []rune) bool {
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

func (s *seq) CharAt(i int) rune {
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
	if idx != nil {
		for _, pos := range idx {
			if pos[1]-pos[0] > end-start {
				end = pos[1]
				start = pos[0]
			}
		}
	}
	return start, end
}

func RandomSequence(alphabet, length int) ([]rune, error) {
	seq := make([]rune, length)
	for i := 0; i < length; i++ {
		switch alphabet {
		case AMINOACIDS:
			seq[i] = stdaminoacid[rand.Intn(len(stdaminoacid))]
		case NUCLEOTIDS:
			seq[i] = stdnucleotides[rand.Intn(len(stdnucleotides))]
		default:
			return nil, errors.New("Unexpected sequence alphabet type")
		}
	}
	return seq, nil
}

// Reverses a sequence
func Reverse(seq []rune) {
	for i, j := 0, len(seq)-1; i < j; i, j = i+1, j-1 {
		seq[i], seq[j] = seq[j], seq[i]
	}
}

// Reverse sequence order
func (s *seq) Reverse() {
	Reverse(s.sequence)
}

// Complement sequence
func Complement(seq []rune) (err error) {

	for i, n := range seq {
		c, ok := complement_nuc_mapping[n]
		if !ok {
			err = fmt.Errorf("Character %c can not be complemented", n)
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
		return fmt.Errorf("Wrong alphabet for Complementing sequence")
	}
	return Complement(s.sequence)
}

func (s *seq) DetectAlphabet() int {
	isaa := true
	isnt := true

	for _, nt := range s.sequence {
		nt = unicode.ToUpper(nt)
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

//NumGaps returns the number of Gaps on the given sequence
func (s *seq) NumGaps() (numgaps int) {
	numgaps = 0
	for _, c := range s.sequence {
		if c == GAP {
			numgaps++
		}
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

// Translates the given sequence start at nucleotide index "phase"
//
// If the sequence is not nucleotidic, then throws an error
// If sequence length is < 3+phase then thropws an error
func (s *seq) Translate(phase int, geneticcode int) (tr Sequence, err error) {
	var buffer bytes.Buffer
	var code map[string]rune

	if code, err = geneticCode(geneticcode); err != nil {
		return
	}

	if err = bufferTranslate(s, phase, code, &buffer); err != nil {
		return
	}

	tr = NewSequence(s.name, []rune(buffer.String()), s.comment)
	return
}

func bufferTranslate(s *seq, phase int, code map[string]rune, buffer *bytes.Buffer) (err error) {
	buffer.Reset()
	if s.DetectAlphabet() != NUCLEOTIDS && s.DetectAlphabet() != BOTH {
		err = fmt.Errorf("Cannot translate this sequence, wrong alphabet")
		return
	}

	if len(s.sequence) < 3+phase {
		err = fmt.Errorf("Cannot translate a sequence with length < 3+phase (%s)", s.name)
		return
	}
	for i := phase; i < len(s.sequence)-2; i += 3 {
		var aa rune = ' '
		var aatmp rune = ' '
		var found bool = false
		// We handle possible IUPAC characters
		codons := GenAllPossibleCodons(s.sequence[i], s.sequence[i+1], s.sequence[i+2])
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
		buffer.WriteRune(aa)
	}
	return
}

func (s *seq) Clone() Sequence {
	seq2 := make([]rune, len(s.sequence))
	copy(seq2, s.sequence)
	return NewSequence(s.name, seq2, s.comment)
}

// GenAllPossibleCodons generates all possible codons given the 3 nucleotides in arguments
// Multiple codons may exist if IUPAC code is employed (R=A|G, etc.).
// The 3 nucleotites in arguments are converted to upper case and U converted to T.
// If one character does not correspond to a known nucleotide in IUPAC code, then
// Returns an empty slice.
//
// For example GenAllPossibleCodons('A','G','N') should return
// {"AGA","AGC","AGG","AGT"}.
func GenAllPossibleCodons(nt1, nt2, nt3 rune) (codons []string) {
	var nt rune
	var codontmp string
	var nts1, nts2, nts3 []rune // possible nts for each nt
	var found bool

	codons = make([]string, 0)
	codonstmp := make([]string, 0)

	nt1 = unicode.ToUpper(nt1)
	nt2 = unicode.ToUpper(nt2)
	nt3 = unicode.ToUpper(nt3)

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
func EqualOrCompatible(nt1, nt2 rune) (ok bool, err error) {
	var ok1, ok2 bool
	var possibilities1, possibilities2 []rune

	nt1 = unicode.ToUpper(nt1)
	nt2 = unicode.ToUpper(nt2)

	possibilities1, ok1 = IupacCode[nt1]
	possibilities2, ok2 = IupacCode[nt2]

	if !ok1 && nt1 != GAP {
		err = fmt.Errorf("Given nucleotide 1 (%c) is not nucleotide", nt1)
		return
	}

	if !ok2 && nt2 != GAP {
		err = fmt.Errorf("Given nucleotide 2 (%c) is not nucleotide", nt2)
		return
	}

	if nt1 == GAP || nt2 == GAP {
		ok = (nt1 == nt2)
		return
	}

	ok = false
	for _, p1 := range possibilities1 {
		for _, p2 := range possibilities2 {
			if p1 == p2 {
				ok = true
				return
			}
		}
	}

	return
}

// NtIUPACDifference returns the cost of the difference between
// the two potentially ambiguous nucleotides.
//
// - if the two nucleotides are identical : returns 0.0
// - if the two nucleotides are different:
//      1) If none are ambigous: returns 1.0
//      2) Otherwise, returns 1-Card(I)/Card(U), I being the
//         intersection of the sets of possible
//         nucleotides of nt1 and nt2, and U being
//         the union of the sets of possible nucleotides
//         of nt1 and nt2.
// For example, if we want to compare Y and S :
// Y = {C | T} and S = {G | C}. Card(I)=1, Card(U)=3, so diff=2/3
//
// Precisions:
// - For N vs. A for example: the difference will be 1-1/4 : 3/4
// - For gaps: Returns diff=1.0
//
// If nt1 or nt2 are not nucleotides, then returns an error
func NtIUPACDifference(nt1, nt2 rune) (diff float64, err error) {
	var possibilities1, possibilities2 []uint8
	var ok1, ok2 bool = true, true
	diff = 0.0

	if nt1 == nt2 {
		return
	}

	diff = 1.0
	nt1 = unicode.ToUpper(nt1)
	nt2 = unicode.ToUpper(nt2)

	possibilities1, ok1 = iupacCodeByte[nt1]

	possibilities2, ok2 = iupacCodeByte[nt2]

	if !ok1 && nt1 != GAP {
		err = fmt.Errorf("Given nucleotide 1 (%c) is not nucleotide", nt1)
		return
	}

	if !ok2 && nt2 != GAP {
		err = fmt.Errorf("Given nucleotide 2 (%c) is not nucleotide", nt2)
		return
	}

	if nt1 == GAP || nt2 == GAP {
		if nt1 == nt2 {
			diff = 0.0
		}
		diff = 1.0
		return
	}

	var leftbyte, rightbyte uint8 = 0, 0
	for _, p1 := range possibilities1 {
		leftbyte |= 1 << p1
	}

	for _, p2 := range possibilities2 {
		rightbyte |= 1 << p2
	}

	inter := bits.OnesCount8(leftbyte & rightbyte)
	union := bits.OnesCount8(leftbyte | rightbyte)
	diff = 1.0 - float64(inter)/float64(union)
	return
}
