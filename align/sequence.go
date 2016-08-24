package align

type Sequence interface {
	Sequence()
	Name()
	Comment()
}

type seq struct {
	name     string // Name of the sequence
	sequence string // Sequence of nucleotides/aa
	comment  string // Comment if any
}

func NewSequence(name string, sequence string, comment string) *seq {
	return &seq{
		name,
		sequence,
		comment,
	}
}

func (s *seq) Sequence() string {
	return s.sequence
}

func (s *seq) Name() string {
	return s.name
}

func (s *seq) Comment() string {
	return s.comment
}
