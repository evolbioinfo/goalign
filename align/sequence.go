package align

import (
	"errors"
	"math/rand"
	"regexp"
	"strings"
)

type Sequence interface {
	Sequence() string
	SequenceChar() []rune
	CharAt(int) rune
	Name() string
	Comment() string
	Length() int
	LongestORF() (start, end int) // Detects the longest ORF in forward strand only
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
	re, _ := regexp.CompilePOSIX("(ATG)(.{3})*(TAA|TGA|TAG)")
	re.Longest()
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
