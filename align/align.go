package align

import (
	"bytes"
	"errors"
	"github.com/fredericlemoine/goalign/io"
	"math/rand"
)

const (
	AMINOACIDS = 0 // Amino acid sequence alphabet
	NUCLEOTIDS = 1 // Nucleotid sequence alphabet

)

type Alignment interface {
	AddSequence(name string, sequence string, comment string) error
	AddSequenceChar(name string, sequence []rune, comment string) error
	GetSequence(name string) (string, bool)
	Iterate(it func(name string, sequence string))
	IterateChar(it func(name string, sequence []rune))
	NbSequences() int
	Length() int
	ShuffleSequences()
	ShuffleSites(rate float64)
	Sample(nb int) (Alignment, error)
	BuildBootstrap() Alignment
	Recombine(rate float64)
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
		io.ExitWithMessage(errors.New("Unexpected sequence alphabet type"))
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
		it(seq.name, string(seq.sequence))
	}
}

func (a *align) IterateChar(it func(name string, sequence []rune)) {
	for _, seq := range a.seqs {
		it(seq.name, seq.sequence)
	}
}

// Adds a sequence to this alignment
func (a *align) AddSequence(name string, sequence string, comment string) error {
	a.AddSequenceChar(name, []rune(sequence), comment)
	return nil
}

func (a *align) AddSequenceChar(name string, sequence []rune, comment string) error {
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

// Just shuffle the order of the sequences in the alignment
// Does not change biological information
func (a *align) ShuffleSequences() {
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

// Shuffles vertically rate sites of the alignment
// randomly
// rate must be >=0 and <=1
func (a *align) ShuffleSites(rate float64) {
	if rate < 0 || rate > 1 {
		io.ExitWithMessage(errors.New("Shuffle site rate must be >=0 and <=1"))
	}
	permutation := rand.Perm(a.Length())
	nb_to_shuffle := int(rate * float64(a.Length()))
	var temp rune
	for i := 0; i < nb_to_shuffle; i++ {
		site := permutation[i]
		var n int = a.NbSequences()
		for n > 1 {
			r := rand.Intn(n)
			n--
			temp = a.seqs[n].sequence[site]
			a.seqs[n].sequence[site] = a.seqs[r].sequence[site]
			a.seqs[r].sequence[site] = temp
		}
	}
}

// recombine a rate of the sequences together
// takes rate/2 seqs and recombine them with the other rate/2 seqs at a random position
// if rate < 0 : does nothing
// if rate > 1 : does nothing
func (a *align) Recombine(rate float64) {
	var nb_to_shuffle, nb_sites int
	var pos int
	var tmpchar rune
	var seq1, seq2 *seq

	if rate < 0 || rate > 1 {
		return
	}
	nb_sites = a.Length()
	nb_to_shuffle = (int)(rate * float64(a.NbSequences()))

	permutation := rand.Perm(a.NbSequences())

	for i := 0; i < int(nb_to_shuffle/2); i++ {
		// We take a random position in the sequences and recombine both
		pos = rand.Intn(nb_sites)
		seq1 = a.seqs[permutation[i]]
		seq2 = a.seqs[permutation[i+(int)(nb_to_shuffle/2)]]
		for pos < nb_sites {
			tmpchar = seq1.sequence[pos]
			seq1.sequence[pos] = seq2.sequence[pos]
			seq2.sequence[pos] = tmpchar
			pos++
		}
	}
}

// Samples randomly a subset of the sequences
// And returns this new alignment
// If nb < 1 or nb > nbsequences returns nil and an error
func (a *align) Sample(nb int) (Alignment, error) {
	if a.NbSequences() < nb || nb < 1 {
		return nil, errors.New("Number of sequences to sample is not compatible with alignment")
	}
	sample := NewAlign(a.alphabet)
	permutation := rand.Perm(a.NbSequences())
	for i := 0; i < nb; i++ {
		seq := a.seqs[permutation[i]]
		sample.AddSequenceChar(seq.name, seq.SequenceChar(), seq.Comment())
	}
	return sample, nil
}

func (a *align) BuildBootstrap() Alignment {
	n := a.Length()
	boot := NewAlign(a.alphabet)
	indices := make([]int, n)
	var buf bytes.Buffer

	for i := 0; i < n; i++ {
		indices[i] = rand.Intn(n)
	}

	for _, seq := range a.seqs {
		buf.Reset()
		for _, indice := range indices {
			buf.WriteRune(seq.sequence[indice])
		}
		boot.AddSequenceChar(seq.name, bytes.Runes(buf.Bytes()), seq.Comment())
	}
	return boot
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
		case 'U', 'O', 'R':
			couldbent = true
		case 'Q', 'E', 'I', 'L', 'F', 'P', 'Z':
			couldbeaa = true
		}
		isaa = isaa && couldbeaa
		isnt = isnt && couldbent
		if !(isaa || isnt) {
			io.ExitWithMessage(errors.New("Unknown character state in alignment : " + string(nt)))
		}
	}

	if isnt {
		return NUCLEOTIDS
	} else {
		return AMINOACIDS
	}
}
