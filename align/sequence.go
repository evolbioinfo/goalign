package align

import (
	"bytes"
	"errors"
	"fmt"
	"math/rand"
	"regexp"
	"strings"
	"unicode"
)

type Sequence interface {
	Sequence() string
	SequenceChar() []rune
	CharAt(int) rune
	Name() string
	SetName(name string)
	Comment() string
	Length() int
	LongestORF() (start, end int) // Detects the longest ORF in forward strand only
	Reverse()
	Complement() error                     // Returns an error if not nucleotide sequence
	Translate(phase int) (Sequence, error) // Translates the sequence using the standard code
	DetectAlphabet() int                   // Try to detect alphabet (nt or aa)
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

// Translates the given sequence start at nucleotide index "phase"
//
// If the sequence is not nucleotidic, then throws an error
// If sequence length is < 3+phase then thropws an error
func (s *seq) Translate(phase int) (tr Sequence, err error) {
	var buffer bytes.Buffer

	if s.DetectAlphabet() != NUCLEOTIDS && s.DetectAlphabet() != BOTH {
		err = fmt.Errorf("Cannot translate this sequence, wrong alphabet")
		return
	}

	if len(s.sequence) < 3+phase {
		err = fmt.Errorf("Cannot translate a sequence with length < 3+phase (%s)", s.name)
		return
	}

	for i := phase; i < len(s.sequence)-2; i += 3 {
		codon := strings.Replace(strings.ToUpper(string(s.sequence[i:i+3])), "U", "T", -1)
		aa, found := standardcode[codon]
		if !found {
			aa = 'X'
		}
		buffer.WriteRune(aa)
	}
	tr = NewSequence(s.name, []rune(buffer.String()), s.comment)
	return
}

func (s *seq) Clone() Sequence {
	seq2 := make([]rune, len(s.sequence))
	copy(seq2, s.sequence)
	return NewSequence(s.name, seq2, s.comment)
}
