package align

import (
	"errors"
	"math/rand"
)

const (
	AMINOACIDS = 0 // Amino acid sequence alphabet
	NUCLEOTIDS = 1 // Nucleotid sequence alphabet

)

type Alignment interface {
	AddSequence(name string, sequence string, comment string) error
	GetSequence(name string) (string, bool)
	Iterate(it func(name string, sequence string))
	NbSequences() int
	Length() int
	Shuffle()
}

type align struct {
	length   int             // Length of alignment
	seqmap   map[string]*seq // Map of sequences
	seqs     []*seq          // Set of sequences (to preserve order)
	alphabet int             // AMINOACIDS OR NUCLEOTIDS
}

func NewAlign(alphabet int) *align {
	switch alphabet {
	case AMINOACIDS, NUCLEOTIDS:
		// OK
	default:
		panic("Unexpected sequence alphabet type")
	}
	return &align{
		-1,
		make(map[string]*seq),
		make([]*seq, 0, 100),
		alphabet,
	}
}

func (a *align) Iterate(it func(name string, sequence string)) {
	for _, seq := range a.seqs {
		it(seq.name, seq.sequence)
	}
}

// Adds a sequence to this alignment
func (a *align) AddSequence(name string, sequence string, comment string) error {
	_, ok := a.seqmap[name]
	if ok {
		return errors.New("Sequence " + name + " already exists in alignment")
	}

	if a.length != -1 && a.length != len(sequence) {
		return errors.New("Sequence " + name + " does not have same length as other sequences")
	}
	a.length = len(sequence)
	seq := NewSequence(name, sequence, comment)
	a.seqmap[name] = seq
	a.seqs = append(a.seqs, seq)
	return nil
}

// If sequence exists in alignment, return true,sequence
// Otherwise, return false,nil
func (a *align) GetSequence(name string) (string, bool) {
	seq, ok := a.seqmap[name]
	return seq.Sequence(), ok
}

// If sequence exists in alignment, return true,sequence
// Otherwise, return false,nil
func (a *align) NbSequences() int {
	return len(a.seqs)
}

func (a *align) Length() int {
	return a.length
}

func (a *align) Shuffle() {
	var temp *seq
	var n int = a.NbSequences()
	for n > 1 {
		r := rand.Intn(n)
		n--
		temp = a.seqs[n]
		a.seqs[n] = a.seqs[r]
		a.seqs[r] = temp
	}
}

func DetectAlphabet(seq string) int {
	isaa := true
	isnt := true

	for _, nt := range seq {
		couldbent := false
		couldbeaa := false
		switch nt {
		case 'A', 'C', 'B', 'G', '?', '-', 'D', 'K', 'S', 'H', 'M', 'N', 'V', 'X', 'T', 'W', 'Y':
			couldbent = true
			couldbeaa = true
		case 'U', 'O':
			couldbent = true
		case 'Q', 'E', 'I', 'L', 'F', 'P', 'Z':
			couldbeaa = true
		}
		isaa = isaa && couldbeaa
		isnt = isnt && couldbent
		if !(isaa || isnt) {
			panic("Unknown character state in alignment : " + string(nt))
		}
	}

	if isnt {
		return NUCLEOTIDS
	} else {
		return AMINOACIDS
	}
}
