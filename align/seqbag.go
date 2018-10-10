package align

import (
	"errors"
	"fmt"
	"log"
	"strings"
	"unicode"

	"github.com/fredericlemoine/goalign/io"
)

type SeqBag interface {
	AddSequence(name string, sequence string, comment string) error
	AddSequenceChar(name string, sequence []rune, comment string) error
	Alphabet() int
	AlphabetStr() string
	AlphabetCharacters() []rune
	AutoAlphabet() // detects and sets alphabet automatically for all the sequences
	GetSequence(name string) (string, bool)
	GetSequenceById(ith int) (string, bool)
	GetSequenceChar(name string) ([]rune, bool)
	GetSequenceCharById(ith int) ([]rune, bool)
	GetSequenceNameById(ith int) (string, bool)
	SetSequenceChar(ithAlign, ithSite int, char rune) error
	Iterate(it func(name string, sequence string))
	IterateChar(it func(name string, sequence []rune))
	IterateAll(it func(name string, sequence []rune, comment string))
	NbSequences() int
}

type seqbag struct {
	seqmap   map[string]*seq // Map of sequences
	seqs     []*seq          // Set of sequences (to preserve order)
	alphabet int             // AMINOACIDS , NUCLEOTIDS or UNKOWN
}

func NewSeqBag(alphabet int) *seqbag {
	switch alphabet {
	case AMINOACIDS, NUCLEOTIDS, UNKNOWN:
		// OK
	default:
		io.ExitWithMessage(errors.New("Unexpected sequence alphabet type"))
	}
	return &seqbag{
		make(map[string]*seq),
		make([]*seq, 0, 100),
		alphabet,
	}
}

// Adds a sequence to this alignment
func (sb *seqbag) AddSequence(name string, sequence string, comment string) error {
	err := sb.AddSequenceChar(name, []rune(sequence), comment)
	return err
}

func (sb *seqbag) AddSequenceChar(name string, sequence []rune, comment string) error {
	_, ok := sb.seqmap[name]
	idx := 0
	tmpname := name
	/* If the sequence name already exists, we add a 4 digit index at the end and print a warning on stderr */
	for ok {
		idx++
		log.Print(fmt.Sprintf("Warning: sequence \"%s\" already exists in alignment, renamed in \"%s_%04d\"", tmpname, name, idx))
		tmpname = fmt.Sprintf("%s_%04d", name, idx)
		_, ok = sb.seqmap[tmpname]
	}
	seq := NewSequence(tmpname, sequence, comment)
	sb.seqmap[tmpname] = seq
	sb.seqs = append(sb.seqs, seq)
	return nil
}

func (sb *seqbag) Alphabet() int {
	return sb.alphabet
}

func (sb *seqbag) AlphabetStr() string {
	switch sb.Alphabet() {
	case NUCLEOTIDS:
		return "nucleotide"
	case AMINOACIDS:
		return "protein"
	default:
		return "unknown"
	}
}

func (sb *seqbag) AlphabetCharacters() (alphabet []rune) {
	if sb.Alphabet() == AMINOACIDS {
		return stdaminoacid
	} else {
		return stdnucleotides
	}
}

// If sequence exists in alignment, return true,sequence
// Otherwise, return false,nil
func (sb *seqbag) GetSequence(name string) (string, bool) {
	seq, ok := sb.seqmap[name]
	if ok {
		return seq.Sequence(), ok
	}
	return "", ok
}

// If sequence exists in alignment, return sequence,true
// Otherwise, return "",false
func (sb *seqbag) GetSequenceById(ith int) (string, bool) {
	if ith >= 0 && ith < sb.NbSequences() {
		return sb.seqs[ith].Sequence(), true
	}
	return "", false
}

// If ith >=0 && i < nbSequences() return sequence,true
// Otherwise, return nil, false
func (sb *seqbag) GetSequenceCharById(ith int) ([]rune, bool) {
	if ith >= 0 && ith < sb.NbSequences() {
		return sb.seqs[ith].SequenceChar(), true
	}
	return nil, false
}

// If sequence exists in alignment, return sequence, true
// Otherwise, return nil,false
func (sb *seqbag) GetSequenceChar(name string) ([]rune, bool) {
	seq, ok := sb.seqmap[name]
	if ok {
		return seq.SequenceChar(), ok
	}
	return nil, false
}

// If ith >=0 && i < nbSequences() return name,true
// Otherwise, return "", false
func (sb *seqbag) GetSequenceNameById(ith int) (string, bool) {
	if ith >= 0 && ith < sb.NbSequences() {
		return sb.seqs[ith].Name(), true
	}
	return "", false
}

// If sequence exists in alignment, return true,sequence
// Otherwise, return false,nil
func (sb *seqbag) NbSequences() int {
	return len(sb.seqs)
}

func (sb *seqbag) SetSequenceChar(ithAlign, ithSite int, char rune) error {
	if ithAlign < 0 || ithAlign >= sb.NbSequences() {
		return errors.New("Sequence index is > number of sequences")
	}
	if ithSite < 0 || ithSite >= sb.seqs[ithAlign].Length() {
		return errors.New("Site index is outside sequence length")
	}

	sb.seqs[ithAlign].sequence[ithSite] = char
	return nil
}

func (sb *seqbag) Iterate(it func(name string, sequence string)) {
	for _, seq := range sb.seqs {
		it(seq.name, string(seq.sequence))
	}
}

func (sb *seqbag) IterateChar(it func(name string, sequence []rune)) {
	for _, seq := range sb.seqs {
		it(seq.name, seq.sequence)
	}
}

func (sb *seqbag) IterateAll(it func(name string, sequence []rune, comment string)) {
	for _, seq := range sb.seqs {
		it(seq.name, seq.sequence, seq.comment)
	}
}

/* It appends the given sequence to the sequence having given name */
func (sb *seqbag) appendToSequence(name string, sequence []rune) error {
	seq, ok := sb.seqmap[name]
	if !ok {
		return errors.New(fmt.Sprintf("Sequence with name %s does not exist in alignment", name))
	}
	seq.sequence = append(seq.sequence, sequence...)
	return nil
}

func (sb *seqbag) AutoAlphabet() {
	isaa := true
	isnt := true

	sb.IterateChar(func(name string, seq []rune) {
		for _, nt := range seq {
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
	})

	if isnt {
		sb.alphabet = NUCLEOTIDS
	} else if isaa {
		sb.alphabet = AMINOACIDS
	} else {
		sb.alphabet = UNKNOWN
	}
}

func DetectAlphabet(seq string) int {
	isaa := true
	isnt := true

	for _, nt := range strings.ToUpper(seq) {
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

		if !(isaa || isnt) {
			return UNKNOWN
		}
	}

	if isnt {
		return NUCLEOTIDS
	} else {
		return AMINOACIDS
	}
}
