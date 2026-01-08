package align

import (
	"bytes"
	"errors"
	"fmt"
	"log"
	"math"
	"math/rand"
	"regexp"
	"sort"
	"strings"
	"unicode"

	"github.com/evolbioinfo/goalign/io"
)

// SeqBag represents a set of unaligned sequences
type SeqBag interface {
	AddSequence(name string, sequence string, comment string) error
	AddSequenceChar(name string, sequence []uint8, comment string) error
	AppendSeqIdentifier(identifier string, right bool)
	Alphabet() int
	SetAlphabet(int) error // Sets the alphabet
	AlphabetStr() string
	AlphabetCharacters() []uint8
	AlphabetCharToIndex(c uint8) int // Returns index of the character (nt or aa) in the AlphabetCharacters() array
	AutoAlphabet()                   // detects and sets alphabet automatically for all the sequences
	DetectAlphabet() (alphabet int)  //  detects the compatible alphabets
	CharStats() map[uint8]int64
	UniqueCharacters() []uint8
	UniqueCharactersWithCase() []uint8
	CharStatsSeq(idx int) (map[uint8]int, error)               // Computes frequency of characters for the given sequence
	CleanNames(namemap map[string]string)                      // Clean sequence names (newick special char)
	Clear()                                                    // Removes all sequences
	CloneSeqBag() (seqs SeqBag, err error)                     // Clones the seqqbag
	Deduplicate(nAsGap bool) (identical [][]string, err error) // Remove duplicate sequences (nAsGap is for considering N/X identical to gaps for sequence comparison)
	FilterLength(minlength, maxlength int) error               // Remove sequences whose length is <minlength or >maxlength
	GetSequence(name string) (string, bool)                    // Get a sequence by names
	GetSequenceById(ith int) (string, bool)
	GetSequenceChar(name string) ([]uint8, bool)
	GetSequenceCharById(ith int) ([]uint8, bool)
	GetSequenceNameById(ith int) (string, bool)
	GetSequenceByName(name string) (Sequence, bool)
	GetSequenceIdByName(name string) (i int) // if the name does not exist, i < 0
	SetSequenceChar(ithAlign, ithSite int, char uint8) error
	// IgnoreIdentical sets the behavior when duplicate names are encountered while building the alignment
	// If ignore is IGNORE_NONE: Does not ignore anything
	// If ignore is IGNORE_NAME: Ignore sequences having the same name (keep the first one whatever their sequence)
	// If ignore is IGNORE_SEQUENCE: Ignore sequences having the same name and the same sequence
	// Otherwise, sets IGNORE_NONE
	IgnoreIdentical(int)
	SampleSeqBag(nb int) (SeqBag, error) // generate a sub sample of the sequences
	Sequence(ith int) (Sequence, bool)
	SequenceByName(name string) (Sequence, bool)
	Identical(SeqBag) bool
	Iterate(it func(name string, sequence string) bool)
	IterateChar(it func(name string, sequence []uint8) bool)
	IterateAll(it func(name string, sequence []uint8, comment string) bool)
	Sequences() []Sequence
	SequencesChan() chan Sequence
	LongestORF(reverse bool) (orf Sequence, err error)
	MaxNameLength() int // maximum sequence name length
	NbSequences() int
	RarefySeqBag(nb int, counts map[string]int) (SeqBag, error) // Take a new rarefied sample taking into accounts weights
	// Removes sequences having >= cutoff gaps, returns number of removed sequences
	RemoveGapSeqs(cutoff float64, ignoreNs bool) int
	// Removes sequences having >= cutoff character, returns number of removed sequences
	RemoveCharacterSeqs(c uint8, cutoff float64, ignoreCase, ignoreGaps, ignoreNs bool) int
	Rename(namemap map[string]string)
	RenameRegexp(regex, replace string, namemap map[string]string) error
	Replace(old, new string, regex bool) error             // Replaces old string with new string in sequences of the alignment
	ReplaceStops(phase int, geneticode int) error          // Replaces stop codons in the given phase using the given genetic code
	ShuffleSequences()                                     // Shuffle sequence order
	String() string                                        // Raw string representation (just write all sequences)
	Translate(phase int, geneticcode int) (err error)      // Translates nt sequence in aa
	ToUpper()                                              // replaces lower case characters by upper case characters
	ToLower()                                              // replaces upper case characters by lower case characters
	ReverseComplement() (err error)                        // Reverse-complements the alignment
	ReverseComplementSequences(name ...string) (err error) // Reverse-complements some sequences in the alignment
	TrimNames(namemap map[string]string, size int) error
	TrimNamesAuto(namemap map[string]string, curid *int) error
	Sort() // Sorts the sequences by name
	Unalign() SeqBag
}

type seqbag struct {
	seqmap          map[string]*seq // Map of sequences
	seqs            []*seq          // Set of sequences (to preserve order)
	ignoreidentical int             // if true, then it won't add the sequence if a sequence with the same name AND same sequence exists
	alphabet        int             // AMINOACIDS , NUCLEOTIDS or UNKOWN
}

func NewSeqBag(alphabet int) *seqbag {
	switch alphabet {
	case AMINOACIDS, NUCLEOTIDS, UNKNOWN:
		// OK
	default:
		io.ExitWithMessage(errors.New("unexpected sequence alphabet type"))
	}
	return &seqbag{
		make(map[string]*seq),
		make([]*seq, 0, 100),
		IGNORE_NONE,
		alphabet,
	}
}

// IgnoreIdentical sets the behavior when duplicates are encountered
// If ignore is IGNORE_NONE: Does not ignore anything
// If ignore is IGNORE_NAME: Ignore sequences having the same name (keep the first one whatever their sequence)
// If ignore is IGNORE_SEQUENCE: Ignore sequences having the same name and the same sequence
// Otherwise, sets IGNORE_NONE
func (sb *seqbag) IgnoreIdentical(ignore int) {
	switch ignore {
	case IGNORE_NONE, IGNORE_NAME, IGNORE_SEQUENCE:
		sb.ignoreidentical = ignore
	default:
		sb.ignoreidentical = IGNORE_NONE
	}
}

// Samples randomly a subset of the sequences
// And returns this new alignment
// If nb < 1 or nb > nbsequences returns nil and an error
func (sb *seqbag) SampleSeqBag(nb int) (sample SeqBag, err error) {
	sample, err = sb.sampleSeqBag(nb)
	return
}

// sampleSeqBag is a private function to allow manipulation of the structure and not the interface
func (sb *seqbag) sampleSeqBag(nb int) (*seqbag, error) {
	if sb.NbSequences() < nb {
		return nil, errors.New("number of sequences to sample is greater than alignment size")
	}
	if nb < 1 {
		return nil, errors.New("cannot sample less than 1 sequence")
	}
	sample := NewSeqBag(sb.alphabet)
	permutation := rand.Perm(sb.NbSequences())
	for i := 0; i < nb; i++ {
		seq := sb.seqs[permutation[i]]
		sample.AddSequenceChar(seq.name, seq.SequenceChar(), seq.Comment())
	}
	return sample, nil
}

// Adds a sequence to this alignment
func (sb *seqbag) AddSequence(name string, sequence string, comment string) error {
	err := sb.AddSequenceChar(name, []uint8(sequence), comment)
	return err
}

// If sb.ignoreidentical is true, then it won't add the sequence if a sequence with the same name AND same sequence
// already exists in the alignment
func (sb *seqbag) AddSequenceChar(name string, sequence []uint8, comment string) error {
	s, ok := sb.seqmap[name]
	idx := 0
	tmpname := name

	// If the sequence name already exists with the same sequence
	// and ignoreidentical is true, then we ignore this sequence
	if ok && sb.ignoreidentical == IGNORE_NAME {
		log.Print(fmt.Sprintf("Warning: sequence name \"%s\" already exists in alignment, ignoring", name))
		return nil
	}

	// If the sequence name already exists with the same sequence
	// and ignoreidentical is true, then we ignore this sequence
	if ok && sb.ignoreidentical == IGNORE_SEQUENCE && s.SameSequence(sequence) {
		log.Print(fmt.Sprintf("Warning: sequence \"%s\" already exists in alignment with the same sequence, ignoring", name))
		return nil
	}

	// Other possibility: we rename the sequence
	// If the sequence name already exists, we add a 4 digit index at the end and print a warning on stderr */
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

// Append a string to all sequence names of the alignment
// If right is true, then append it to the right of each names,
// otherwise, appends it to the left
func (sb *seqbag) AppendSeqIdentifier(identifier string, right bool) {
	if len(identifier) != 0 {
		for _, seq := range sb.seqs {
			if right {
				seq.name = seq.name + identifier
			} else {
				seq.name = identifier + seq.name
			}
		}
	}
}

func (sb *seqbag) Alphabet() int {
	return sb.alphabet
}

// SetAlphabet sets the alphabet of the current sequences
// alphabet can be align.AMINOACIDS or align.NUCLEOTIDS
// otherwise returns an error
func (sb *seqbag) SetAlphabet(alphabet int) (err error) {
	alpha := sb.DetectAlphabet()
	if alpha == UNKNOWN {
		err = fmt.Errorf("input sequence alphabet unknown")
		return
	}

	if alphabet == NUCLEOTIDS {
		if alpha == NUCLEOTIDS || alpha == BOTH {
			sb.alphabet = NUCLEOTIDS
		} else {
			err = fmt.Errorf("given alphabet is not compatible with input sequences")
			return
		}
	} else if alphabet == AMINOACIDS {
		if alpha == AMINOACIDS || alpha == BOTH {
			sb.alphabet = AMINOACIDS
		} else {
			err = fmt.Errorf("given alphabet is not compatible with input sequences")
			return
		}
	} else {
		err = fmt.Errorf("given alphabet can not be used in alignment")
		return
	}

	return
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

func (sb *seqbag) AlphabetCharacters() (alphabet []uint8) {
	if sb.Alphabet() == AMINOACIDS {
		return stdaminoacid
	} else {
		return stdnucleotides
	}
}

// Returns index of the character (nt or aa) in the AlphabetCharacters() array
// If character does not exist or alphabet is unkown, then returns -1
func (sb *seqbag) AlphabetCharToIndex(c uint8) int {
	switch sb.Alphabet() {
	case AMINOACIDS:
		if c, err := AA2Index(c); err != nil {
			return -1
		} else {
			return c
		}
	case NUCLEOTIDS:
		if c, err := Nt2Index(c); err != nil {
			return -1
		} else {
			return c
		}
	default:
		return -1
	}
}

// Removes spaces and tabs at beginning and end of sequence names
// and replaces newick special characters \s\t()[];,.: by "-"
func (sb *seqbag) CleanNames(namemap map[string]string) {
	firstlast := regexp.MustCompile(`(^[\s\t]+|[\s\t]+$)`)
	inside := regexp.MustCompile(`[\|\s\t,\[\]\(\),;\.:]+`)

	for _, seq := range sb.seqs {
		old := seq.name
		seq.name = firstlast.ReplaceAllString(seq.name, "")
		seq.name = inside.ReplaceAllString(seq.name, "-")
		if namemap != nil {
			namemap[old] = seq.name
		}
	}
}

// Removes all the sequences from the seqbag
func (sb *seqbag) Clear() {
	sb.seqmap = make(map[string]*seq)
	sb.seqs = make([]*seq, 0, 100)
}

func (sb *seqbag) CloneSeqBag() (SeqBag, error) {
	c := NewSeqBag(sb.Alphabet())
	c.IgnoreIdentical(sb.ignoreidentical)
	var err error
	sb.IterateAll(func(name string, sequence []uint8, comment string) bool {
		newseq := make([]uint8, 0, len(sequence))
		newseq = append(newseq, sequence...)
		err = c.AddSequenceChar(name, newseq, comment)
		return err != nil
	})
	return c, err
}

// This function removes sequences that are duplicates of other
// It keeps one copy of each sequence, with the name of the first
// found.
// nAsGap: if true, then considers N characters / X characters as identical to GAPs for sequence comparison
//
// As output, identical contains a slice of identical sequence names
// ex: identical[0] is a slice of identical sequence names
//
// It modifies input alignment.
func (sb *seqbag) Deduplicate(nAsGap bool) (identical [][]string, err error) {
	var compareString string
	oldseqs := sb.seqs
	sb.Clear()
	identical = make([][]string, 0)

	// Stores the index of the identical sequence group
	// of a given sequence in the "identical" slice of slice
	seqs := make(map[string]int)
	for _, seq := range oldseqs {
		s := string(seq.sequence)

		// We create a temp sequence with N/X replaced by GAP
		if nAsGap {
			if sb.Alphabet() == AMINOACIDS {
				compareString = strings.ReplaceAll(s, string(ALL_AMINO), string(GAP))
			} else if sb.Alphabet() == NUCLEOTIDS {
				compareString = strings.ReplaceAll(s, string(ALL_NUCLE), string(GAP))
			} else {
				compareString = s
			}
		} else {
			compareString = s
		}
		// If the group does not exist
		if i, ok := seqs[compareString]; !ok {
			if err = sb.AddSequence(seq.name, s, seq.comment); err != nil {
				return
			}
			identical = append(identical, []string{seq.name})
			seqs[compareString] = len(identical) - 1
		} else {
			identical[i] = append(identical[i], seq.name)
		}
	}
	return
}

// FilterLength removes sequences whose length is <minlength or >maxlength
// If maxlength < 0 : only considers minlength
// If minlength < 0 : only considers maxlength
//
// The original alignment is directly modified. To work on a copy of the original
// alignment, think about SeqBag.CloneSeqBag()
func (sb *seqbag) FilterLength(minlength, maxlength int) (err error) {
	oldseqs := sb.seqs
	sb.Clear()
	for _, seq := range oldseqs {
		if (minlength >= 0 && seq.Length() >= minlength) || (maxlength > 0 && seq.Length() <= maxlength) {
			if err = sb.AddSequenceChar(seq.name, seq.sequence, seq.comment); err != nil {
				return
			}
		}
	}
	return
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

// If sequence exists in alignment, return true,sequence
// Otherwise, return false,nil
func (sb *seqbag) GetSequenceByName(name string) (Sequence, bool) {
	if seq, ok := sb.seqmap[name]; ok {
		return seq, true
	}
	return nil, false
}

// If sequence exists in alignment, return its index
// Otherwise, return -1
func (sb *seqbag) GetSequenceIdByName(name string) int {
	var s *seq
	var i int
	for i, s = range sb.seqs {
		if s.name == name {
			return i
		}
	}
	return -1
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
func (sb *seqbag) GetSequenceCharById(ith int) ([]uint8, bool) {
	if ith >= 0 && ith < sb.NbSequences() {
		return sb.seqs[ith].SequenceChar(), true
	}
	return nil, false
}

// If sequence exists in alignment, return sequence, true
// Otherwise, return nil,false
func (sb *seqbag) GetSequenceChar(name string) ([]uint8, bool) {
	seq, ok := sb.seqmap[name]
	if ok {
		return seq.SequenceChar(), ok
	}
	return nil, false
}

// If sequence exists in alignment, return sequence,true
// Otherwise, return "",false
func (sb *seqbag) Sequence(ith int) (Sequence, bool) {
	if ith >= 0 && ith < sb.NbSequences() {
		return sb.seqs[ith], true
	}
	return nil, false
}

// If sequence exists in alignment, return sequence,true
// Otherwise, return "",false
func (sb *seqbag) SequenceByName(name string) (Sequence, bool) {
	seq, ok := sb.seqmap[name]
	if ok {
		return seq, ok
	}
	return nil, ok
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

func (sb *seqbag) SetSequenceChar(ithAlign, ithSite int, char uint8) error {
	if ithAlign < 0 || ithAlign >= sb.NbSequences() {
		return errors.New("sequence index is > number of sequences")
	}
	if ithSite < 0 || ithSite >= sb.seqs[ithAlign].Length() {
		return errors.New("site index is outside sequence length")
	}

	sb.seqs[ithAlign].sequence[ithSite] = char
	return nil
}

// Returns true if:
//
//   - sb and comp have the same number of sequences &&
//   - each sequence in sb have a sequence in comp with the same name
//     and the same sequence
//
// Identical seqbags may have sequences in different order
func (sb *seqbag) Identical(comp SeqBag) bool {
	if sb.NbSequences() != comp.NbSequences() {
		return false
	}

	for _, seq := range sb.seqs {
		seq2, ok := comp.GetSequence(seq.name)
		if !ok {
			return false
		}
		if string(seq.sequence) != seq2 {
			return false
		}
	}

	return true
}

func (sb *seqbag) Iterate(it func(name string, sequence string) bool) {
	var stop bool = false
	for _, seq := range sb.seqs {
		if stop = it(seq.name, string(seq.sequence)); stop {
			return
		}
	}
}

func (sb *seqbag) IterateChar(it func(name string, sequence []uint8) bool) {
	var stop bool = false
	for _, seq := range sb.seqs {
		if stop = it(seq.name, seq.sequence); stop {
			return
		}
	}
}

func (sb *seqbag) IterateAll(it func(name string, sequence []uint8, comment string) bool) {
	var stop bool = false
	for _, seq := range sb.seqs {
		if stop = it(seq.name, seq.sequence, seq.comment); stop {
			return
		}
	}
}

func (sb *seqbag) Sequences() (seqs []Sequence) {
	seqs = make([]Sequence, len(sb.seqs))
	for i, s := range sb.seqs {
		seqs[i] = s
	}
	return seqs
}

func (sb *seqbag) SequencesChan() (seqs chan Sequence) {
	seqs = make(chan Sequence, 50)
	go func() {
		for _, s := range sb.seqs {
			seqs <- s
		}
		close(seqs)
	}()
	return
}

/* It appends the given sequence to the sequence having given name */
func (sb *seqbag) appendToSequence(name string, sequence []uint8) error {
	seq, ok := sb.seqmap[name]
	if !ok {
		return fmt.Errorf("Sequence with name %s does not exist in alignment", name)
	}
	seq.sequence = append(seq.sequence, sequence...)
	return nil
}

// Detects the alphabets compatible with the alignment
// can be align.BOTH, align.NUCLEOTIDS, align.AMINOACIDS, or align.UNKNOWN
func (sb *seqbag) DetectAlphabet() (alphabet int) {
	isaa := true
	isnt := true

	sb.IterateChar(func(name string, seq []uint8) bool {
		for _, nt := range seq {
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
		return false
	})

	if isnt {
		if isaa {
			alphabet = BOTH
		} else {
			alphabet = NUCLEOTIDS
		}
	} else {
		if isaa {
			alphabet = AMINOACIDS
		} else {
			alphabet = UNKNOWN
		}
	}
	return
}

func (sb *seqbag) AutoAlphabet() {
	alphabet := sb.DetectAlphabet()

	if alphabet == BOTH || alphabet == NUCLEOTIDS {
		sb.alphabet = NUCLEOTIDS
	} else if alphabet == AMINOACIDS {
		sb.alphabet = AMINOACIDS
	} else {
		sb.alphabet = UNKNOWN
	}
}

// Returns the distribution of all characters
func (sb *seqbag) CharStats() (chars map[uint8]int64) {
	chars = make(map[uint8]int64)
	present := make([]int, 130)

	for _, seq := range sb.seqs {
		for _, r := range seq.sequence {
			present[unicode.ToUpper(rune(r))]++
		}
	}

	for i, r := range present {
		if r > 0 {
			chars[uint8(i)] = int64(r)
		}
	}

	return
}

// Returns the distribution of all characters (all uppercase)
func (sb *seqbag) UniqueCharacters() (chars []uint8) {
	chars = make([]uint8, 0, 40)
	present := make([]bool, 130)

	for _, seq := range sb.seqs {
		for _, r := range seq.sequence {
			present[unicode.ToUpper(rune(r))] = true
		}
	}

	for i, r := range present {
		if r {
			chars = append(chars, uint8(i))
		}
	}

	return
}

// Returns the distribution of all characters keeping the original case
func (sb *seqbag) UniqueCharactersWithCase() (chars []uint8) {
	chars = make([]uint8, 0, 40)
	present := make([]bool, 130)

	for _, seq := range sb.seqs {
		for _, r := range seq.sequence {
			present[rune(r)] = true
		}
	}

	for i, r := range present {
		if r {
			chars = append(chars, uint8(i))
		}
	}

	return
}

// CharStatsSeq Returns the frequency of all characters in the sequence
// identified by the given index.
// If the sequence with the given index does not exist, then returns an error
func (sb *seqbag) CharStatsSeq(idx int) (outmap map[uint8]int, err error) {
	var seq []uint8
	var ok bool

	outmap = make(map[uint8]int)

	if seq, ok = sb.GetSequenceCharById(idx); !ok {
		err = fmt.Errorf("Sequence with id %d does not exist", idx)
	} else {
		for _, r := range seq {
			outmap[uint8(unicode.ToUpper(rune(r)))]++
		}
	}

	return
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

	if isnt && isaa {
		return BOTH
	} else if isnt {
		return NUCLEOTIDS
	} else {
		return AMINOACIDS
	}
}

/*
Each sequence in the alignment has an associated number of occurence. The sum s of the counts
represents the number of sequences in the underlying initial dataset.

The goal is to downsample (rarefy) the initial dataset, by sampling n sequences
from s (n<s), and taking the alignment corresponding to this new sample, i.e by
taking only unique (different) sequences from it.

Parameters are:
* nb: the number of sequences to sample from the underlying full dataset (different
from the number of sequences in the output alignment)
* counts: counts associated to each sequence (if the count of a sequence is missing, it
is considered as 0, if the count of an unkown sequence is present, it will return an error).

	Sum of counts of all sequences must be > n.
*/
func (sb *seqbag) RarefySeqBag(nb int, counts map[string]int) (sample SeqBag, err error) {
	sample, err = sb.rarefySeqBag(nb, counts)
	return
}

// rarefySeqBag is a private function to allow manipulation of the structure and not the interface
func (sb *seqbag) rarefySeqBag(nb int, counts map[string]int) (sample *seqbag, err error) {
	// Sequences that will be selected
	selected := make(map[string]bool)
	total := 0
	// We copy the count map to modify it
	tmpcounts := make(map[string]int)
	tmpcountskeys := make([]string, len(counts))
	i := 0
	for k, v := range counts {
		tmpcountskeys[i] = k
		if v <= 0 {
			err = fmt.Errorf("Sequence counts must be positive")
			return
		}
		if _, ok := sb.GetSequenceChar(k); !ok {
			err = fmt.Errorf("Sequence %s does not exist in the alignment", k)
			return
		}
		tmpcounts[k] = v
		total += v
		i++
	}

	sort.Strings(tmpcountskeys)

	if nb >= total {
		err = fmt.Errorf("number of sequences to sample %d is >= sum of the counts %d", nb, total)
		return
	}

	// We sample a new sequence nb times
	for i := 0; i < nb; i++ {
		proba := 0.0
		// random num
		unif := rand.Float64()
		for idk, k := range tmpcountskeys {
			v, ok := tmpcounts[k]
			if !ok {
				err = fmt.Errorf("no sequence named %s is present in the tmp count map", k)
				return
			}
			proba += float64(v) / float64(total)
			if unif < proba {
				selected[k] = true
				if v-1 == 0 {
					delete(tmpcounts, k)
					tmpcountskeys = append(tmpcountskeys[:idk], tmpcountskeys[idk+1:]...)
				} else {
					tmpcounts[k] = v - 1
				}
				break
			}
		}
		total--
	}

	sample = NewSeqBag(sb.alphabet)
	sb.IterateAll(func(name string, sequence []uint8, comment string) bool {
		if _, ok := selected[name]; ok {
			sample.AddSequenceChar(name, sequence, comment)
		}
		return false
	})

	return
}

// Removes sequences constituted of [cutoff*100%,100%] Gaps
// Exception fo a cutoff of 0: does not remove sequences with 0% gaps
// Cutoff must be between 0 and 1, otherwise set to 0.
// 0 means that sequences with > 0 gaps will be removed
// other cutoffs : ]0,1] mean that sequences with >= cutoff gaps will be removed
//
// Returns the number of removed sequences
func (sb *seqbag) RemoveGapSeqs(cutoff float64, ignoreNs bool) int {
	return sb.RemoveCharacterSeqs(GAP, cutoff, false, false, ignoreNs)
}

// RemoveCharacterSeqs Removes sequences constituted of [cutoff*100%,100%] of the given character
// Exception fo a cutoff of 0: does not remove sequences with 0% of this character
// Cutoff must be between 0 and 1, otherwise set to 0.
// 0 means that sequences with > 0 of the given character will be removed
// other cutoffs : ]0,1] mean that positions with >= cutoff of this character will be removed
//
// if ignoreCase then the search is case insensitive
// if ignoreGaps is true, then gaps are not taken into account
// if ignoreNs is true, then Ns are not taken into account
//
// Returns the number of removed sequences
func (sb *seqbag) RemoveCharacterSeqs(c uint8, cutoff float64, ignoreCase, ignoreGaps, ignoreNs bool) int {
	var nbseqs int
	var total int
	if cutoff < 0 || cutoff > 1 {
		cutoff = 0
	}
	oldseqs := sb.seqs
	sb.Clear()
	nbremoved := 0

	all := uint8(ALL_NUCLE)
	if sb.Alphabet() == AMINOACIDS {
		all = uint8(ALL_AMINO)
	}
	allc := uint8(unicode.ToLower(rune(all)))

	for _, seq := range oldseqs {
		nbseqs = 0
		total = 0

		for site := 0; site < seq.Length(); site++ {
			if (seq.sequence[site] == c) || (ignoreCase && unicode.ToLower(rune(seq.sequence[site])) == unicode.ToLower(rune(c))) {
				nbseqs++
			}
			// If we exclude gaps and it is a gap: we do nothing
			// or if we exclude Ns and it is a N: we do nothing
			if !(ignoreGaps && seq.sequence[site] == uint8(GAP)) && !(ignoreNs && (seq.sequence[site] == all || seq.sequence[site] == allc)) {
				total++
			}
		}
		if (cutoff > 0.0 && float64(nbseqs) >= cutoff*float64(total)) || (cutoff == 0 && nbseqs > 0) {
			nbremoved++
		} else {
			sb.AddSequenceChar(seq.name, seq.sequence, seq.comment)
		}
	}

	return nbremoved
}

// This function renames sequences of the alignment based on the map in argument
// If a name in the map does not exist in the alignment, does nothing
// If a sequence in the alignment does not have a name in the map: does nothing
func (sb *seqbag) Rename(namemap map[string]string) {
	for seq := 0; seq < sb.NbSequences(); seq++ {
		newname, ok := namemap[sb.seqs[seq].name]
		if ok {
			sb.seqs[seq].name = newname
		}
		// else {
		// 	io.PrintMessage("Sequence " + a.seqs[seq].name + " not present in the map file")
		// }
	}
}

// Shuffle the order of the sequences in the alignment
// Does not change biological information
func (sb *seqbag) ShuffleSequences() {
	var temp *seq
	var n int = sb.NbSequences()
	for n > 1 {
		r := rand.Intn(n)
		n--
		temp = sb.seqs[n]
		sb.seqs[n] = sb.seqs[r]
		sb.seqs[r] = temp
	}
}

// This function renames sequences of the alignment based on the given regex and replace strings
func (sb *seqbag) RenameRegexp(regex, replace string, namemap map[string]string) error {
	r, err := regexp.Compile(regex)
	if err != nil {
		return err
	}
	for seq := 0; seq < sb.NbSequences(); seq++ {
		newname := r.ReplaceAllString(sb.seqs[seq].name, replace)
		namemap[sb.seqs[seq].name] = newname
		sb.seqs[seq].name = newname
	}
	return nil
}

// Replace an old string in sequences by a new string
// It may be a regexp
//
// - If the regex is malformed, returns an error
func (sb *seqbag) Replace(old, new string, regex bool) (err error) {
	var r *regexp.Regexp

	if regex {
		r, err = regexp.Compile(old)
		if err != nil {
			return err
		}
		for seq := 0; seq < sb.NbSequences(); seq++ {
			newseq := []uint8(r.ReplaceAllString(string(sb.seqs[seq].sequence), new))
			sb.seqs[seq].sequence = newseq
		}
	} else {
		for seq := 0; seq < sb.NbSequences(); seq++ {
			newseq := strings.Replace(string(sb.seqs[seq].sequence), old, new, -1)
			sb.seqs[seq].sequence = []uint8(newseq)
		}
	}
	return nil
}

// Sorts sequences by name
func (sb *seqbag) Sort() {
	names := make([]string, len(sb.seqs))

	// Get sequence names
	for i, seq := range sb.seqs {
		names[i] = seq.Name()
	}

	// Sort names
	sort.Strings(names)
	for i, n := range names {
		s := sb.seqmap[n]
		sb.seqs[i] = s
	}
}

// Replace an old string in sequences by a new string
// It may be a regexp
// Uses the given genetic code
func (sb *seqbag) ReplaceStops(phase int, geneticcode int) (err error) {
	var code map[string]uint8
	var aa uint8

	if sb.Alphabet() != NUCLEOTIDS {
		err = errors.New("wrong alphabet, cannot replace stop codons")
		return
	}

	if code, err = geneticCode(geneticcode); err != nil {
		return
	}

	for seq := 0; seq < sb.NbSequences(); seq++ {
		s := sb.seqs[seq]
		for codon := phase; codon < sb.seqs[seq].Length()-5; codon += 3 {
			aa = translateCodon(s.sequence[codon], s.sequence[codon+1], s.sequence[codon+2], code)
			if aa == OTHER {
				s.sequence[codon] = ALL_NUCLE
				s.sequence[codon+1] = ALL_NUCLE
				s.sequence[codon+2] = ALL_NUCLE
			}
		}
	}
	return nil
}

/*
Translates the given sequences in aa, in the given phase (0,1,2). All sequences
are consifered being in the same phase.

- if the phase is 1 or 2 : it will remove the first or the 2 first characters
- if the phase is -1: it will translate in the 3 phases, and append the suffix _<phase> to all
sequence names. At the end 3x more sequences in the seqbag.
- if the alphabet is not NUCLEOTIDES: returns an error

The seqbag is cleared and old sequences are replaced with aminoacid sequences
*/
func (sb *seqbag) Translate(phase int, geneticcode int) (err error) {
	var oldseqs []*seq
	var buffer bytes.Buffer
	var firststart, laststart int
	var name string
	var suffix bool
	var code map[string]uint8

	if code, err = geneticCode(geneticcode); err != nil {
		return
	}

	if sb.Alphabet() != NUCLEOTIDS {
		err = errors.New("wrong alphabet, cannot translate to AA")
		return
	}

	oldseqs = sb.seqs
	sb.Clear()

	firststart = phase
	laststart = phase
	suffix = false
	if phase == -1 {
		firststart = 0
		laststart = 2
		suffix = true
	}

	for _, seq := range oldseqs {
		name = seq.name
		// We may translate in several phases (if phase==-1)
		for phase = firststart; phase <= laststart; phase++ {
			if suffix {
				name = fmt.Sprintf("%s_%d", seq.name, phase)
			}

			if err = bufferTranslate(seq, phase, code, &buffer); err != nil {
				return
			}

			if err = sb.AddSequence(name, buffer.String(), seq.comment); err != nil {
				return
			}
		}
	}
	sb.AutoAlphabet()

	return
}

/*
ToUpper replaces all lowercase characters by upper case characters in the input sequences.

The seqbag is cleared and old sequences are replaced with upper case sequences
*/
func (sb *seqbag) ToUpper() {
	for _, seq := range sb.seqs {
		for i, c := range seq.sequence {
			seq.sequence[i] = uint8(unicode.ToUpper(rune(c)))
		}
	}
}

/*
ToLower replaces all uppercase characters by lower case characters in the input sequences.

The seqbag is cleared and old sequences are replaced with lower case sequences
*/
func (sb *seqbag) ToLower() {
	for _, seq := range sb.seqs {
		for i, c := range seq.sequence {
			seq.sequence[i] = uint8(unicode.ToLower(rune(c)))
		}
	}
}

/*
ReverseComplement reverse complements all input seuqences.
- if the alphabet is not NUCLEOTIDES: returns an error
- IUPAC characters are supported
*/
func (sb *seqbag) ReverseComplement() (err error) {
	if sb.Alphabet() != NUCLEOTIDS {
		err = errors.New("wrong alphabet, cannot reverse complement")
		return
	}

	for _, seq := range sb.seqs {
		if err = Complement(seq.sequence); err != nil {
			return
		}
		Reverse(seq.sequence)
	}

	return
}

/*
ReverseComplement reverse complements all input seuqences.
- if the alphabet is not NUCLEOTIDES: returns an error
- IUPAC characters are supported
*/
func (sb *seqbag) ReverseComplementSequences(names ...string) (err error) {
	if sb.Alphabet() != NUCLEOTIDS {
		err = errors.New("wrong alphabet, cannot reverse complement")
		return
	}

	for _, name := range names {
		s, found := sb.SequenceByName(name)
		if found {
			if err = Complement(s.SequenceChar()); err != nil {
				return
			}
			Reverse(s.SequenceChar())
		}
	}

	return
}

// Translate sequences in 3 phases (or 6 phases if reverse strand is true)
// And return the longest orf found
func (sb *seqbag) LongestORF(reverse bool) (orf Sequence, err error) {
	var beststart, bestend int
	var start, end int
	var bestseq Sequence
	var found bool
	var name string

	found = false
	// Search for the longest orf in all sequences
	for _, seq := range sb.seqs {
		start, end = seq.LongestORF()
		if start != -1 && end-start > bestend-beststart {
			beststart, bestend = start, end
			bestseq = seq
			name = seq.name
			found = true
		}
		if reverse {
			rev := seq.Clone()
			rev.Reverse()
			rev.Complement()
			start, end = rev.LongestORF()
			if start != -1 && end-start > bestend-beststart {
				beststart, bestend = start, end
				bestseq = rev
				name = seq.name
				found = true
			}
		}
	}

	if !found {
		err = fmt.Errorf("no ORF has been found on any sequence")
		return
	}

	// log.Print("Longest ORF found in sequence ", bestseq.Name())
	// log.Print(string(bestseq.SequenceChar()[beststart:bestend]))
	orf = NewSequence(name, bestseq.SequenceChar()[beststart:bestend], "")
	return
}

func (sb *seqbag) MaxNameLength() (max int) {
	max = 0
	for _, s := range sb.seqs {
		if len(s.Name()) > max {
			max = len(s.Name())
		}
	}
	return
}

// Shorten sequence names to the given size. If duplicates are generated
// then add an identifier.
//
// The map in argument is updated with new oldname=>newname key values
func (sb *seqbag) TrimNames(namemap map[string]string, size int) error {
	shortmap := make(map[string]bool)
	if math.Pow10(size-2) < float64(sb.NbSequences()) {
		return fmt.Errorf("new name size (%d) does not allow to identify that amount of sequences (%d)",
			size-2, sb.NbSequences())
	}
	// If previous short names, we take them into account for uniqueness
	for _, v := range namemap {
		shortmap[v] = true
	}
	for _, seq := range sb.seqs {
		newname, ok := namemap[seq.Name()]
		if !ok {
			newname = seq.Name()
			newname = strings.Replace(newname, ":", "", -1)
			newname = strings.Replace(newname, "_", "", -1)
			if len(newname) >= size-2 {
				newname = newname[0 : size-2]
			} else {
				for m := 0; m < (size - 2 - len(newname)); m++ {
					newname = newname + "x"
				}
			}
			id := 1
			_, ok2 := shortmap[fmt.Sprintf("%s%02d", newname, id)]
			for ok2 {
				id++
				if id > 99 {
					return errors.New("More than 100 identical short names (" + newname + "), cannot shorten the names")
				}
				_, ok2 = shortmap[fmt.Sprintf("%s%02d", newname, id)]
			}
			newname = fmt.Sprintf("%s%02d", newname, id)
			shortmap[newname] = true
			namemap[seq.Name()] = newname
		}
		delete(sb.seqmap, seq.name)
		seq.name = newname
		sb.seqmap[seq.name] = seq
	}

	return nil
}

// namemap : map that is updated during the process
// curid : pointer to the current identifier to use for next seq names
// it is incremented in the function.
func (sb *seqbag) TrimNamesAuto(namemap map[string]string, curid *int) (err error) {
	length := int(math.Ceil(math.Log10(float64(sb.NbSequences() + 1))))
	for _, seq := range sb.seqs {
		newname, ok := namemap[seq.Name()]
		if !ok {
			newname = fmt.Sprintf(fmt.Sprintf("S%%0%dd", (length)), *curid)
			namemap[seq.Name()] = newname
			(*curid)++
			// In case of several alignments to rename,
			// The number of necessary digits may be updated
			length = int(math.Ceil(math.Log10(float64(*curid + 1))))
		}
		seq.name = newname
	}
	return
}

func (sb *seqbag) Unalign() (unal SeqBag) {
	unal = NewSeqBag(sb.Alphabet())

	for _, seq := range sb.seqs {
		unal.AddSequence(seq.name, strings.Replace(string(seq.sequence), "-", "", -1), seq.comment)
	}
	return
}

func (sb *seqbag) String() string {
	var buffer bytes.Buffer
	buffer.WriteString("\n")
	for _, seq := range sb.seqs {
		buffer.WriteString(seq.name + ":")
		buffer.WriteString(string(seq.sequence))
		buffer.WriteRune('\n')
	}
	return buffer.String()
}
