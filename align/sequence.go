package align

import (
	"errors"
	"math/rand"
)

type Sequence interface {
	Sequence() string
	SequenceChar() []rune
	Name() string
	Comment() string
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

func (s *seq) Name() string {
	return s.name
}

func (s *seq) Comment() string {
	return s.comment
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
