package align

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
