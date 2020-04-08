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
	AddSequenceChar(name string, sequence []rune, comment string) error
	AppendSeqIdentifier(identifier string, right bool)
	Alphabet() int
	AlphabetStr() string
	AlphabetCharacters() []rune
	AlphabetCharToIndex(c rune) int // Returns index of the character (nt or aa) in the AlphabetCharacters() array
	AutoAlphabet()                  // detects and sets alphabet automatically for all the sequences
	CharStats() map[rune]int64
	CleanNames(namemap map[string]string)           // Clean sequence names (newick special char)
	Clear()                                         // Removes all sequences
	CloneSeqBag() (seqs SeqBag, err error)          // Clones the seqqbag
	Deduplicate() (identical [][]string, err error) // Remove duplicate sequences
	FilterLength(minlength, maxlength int) error    // Remove sequences whose length is <minlength or >maxlength
	GetSequence(name string) (string, bool)         // Get a sequence by names
	GetSequenceById(ith int) (string, bool)
	GetSequenceChar(name string) ([]rune, bool)
	GetSequenceCharById(ith int) ([]rune, bool)
	GetSequenceNameById(ith int) (string, bool)
	SetSequenceChar(ithAlign, ithSite int, char rune) error
	IgnoreIdentical(bool)                // if true, then it won't add the sequence if a sequence with the same name AND same sequence exists
	SampleSeqBag(nb int) (SeqBag, error) // generate a sub sample of the sequences
	Sequence(ith int) (Sequence, bool)
	SequenceByName(name string) (Sequence, bool)
	Identical(SeqBag) bool
	Iterate(it func(name string, sequence string) bool)
	IterateChar(it func(name string, sequence []rune) bool)
	IterateAll(it func(name string, sequence []rune, comment string) bool)
	Sequences() []Sequence
	SequencesChan() chan Sequence
	LongestORF(reverse bool) (orf Sequence, err error)
	MaxNameLength() int // maximum sequence name length
	NbSequences() int
	RarefySeqBag(nb int, counts map[string]int) (SeqBag, error) // Take a new rarefied sample taking into accounts weights
	Rename(namemap map[string]string)
	RenameRegexp(regex, replace string, namemap map[string]string) error
	Replace(old, new string, regex bool) error        // Replaces old string with new string in sequences of the alignment
	ShuffleSequences()                                // Shuffle sequence order
	String() string                                   // Raw string representation (just write all sequences)
	Translate(phase int, geneticcode int) (err error) // Translates nt sequence in aa
	TrimNames(namemap map[string]string, size int) error
	TrimNamesAuto(namemap map[string]string, curid *int) error
	Sort() // Sorts the sequences by name
	Unalign() SeqBag
}

type seqbag struct {
	seqmap          map[string]*seq // Map of sequences
	seqs            []*seq          // Set of sequences (to preserve order)
	ignoreidentical bool            // if true, then it won't add the sequence if a sequence with the same name AND same sequence exists
	alphabet        int             // AMINOACIDS , NUCLEOTIDS or UNKOWN
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
		false,
		alphabet,
	}
}

func (sb *seqbag) IgnoreIdentical(ignoreidentical bool) {
	sb.ignoreidentical = ignoreidentical
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
		return nil, errors.New("Number of sequences to sample is greater than alignment size")
	}
	if nb < 1 {
		return nil, errors.New("Cannot sample less than 1 sequence")
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
	err := sb.AddSequenceChar(name, []rune(sequence), comment)
	return err
}

// If sb.ignoreidentical is true, then it won't add the sequence if a sequence with the same name AND same sequence
// already exists in the alignment
func (sb *seqbag) AddSequenceChar(name string, sequence []rune, comment string) error {
	s, ok := sb.seqmap[name]
	idx := 0
	tmpname := name

	// If the sequence name already exists with the same sequence
	// and ignoreidentical is true, then we ignore this sequence
	if ok && sb.ignoreidentical && s.SameSequence(sequence) {
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

// Returns index of the character (nt or aa) in the AlphabetCharacters() array
// If character does not exist or alphabet is unkown, then returns -1
func (sb *seqbag) AlphabetCharToIndex(c rune) int {
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
	firstlast := regexp.MustCompile("(^[\\s\\t]+|[\\s\\t]+$)")
	inside := regexp.MustCompile("[\\|\\s\\t,\\[\\]\\(\\),;\\.:]+")

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
	sb.IterateAll(func(name string, sequence []rune, comment string) bool {
		newseq := make([]rune, 0, len(sequence))
		newseq = append(newseq, sequence...)
		err = c.AddSequenceChar(name, newseq, comment)
		return err != nil
	})
	return c, err
}

// This function removes sequences that are duplicates of other
// It keeps one copy of each sequence, with the name of the first
// found.
//
// As output, identical contains a slice of identical sequence names
// ex: identical[0] is a slice of identical sequence names
//
// It modifies input alignment.
func (sb *seqbag) Deduplicate() (identical [][]string, err error) {
	oldseqs := sb.seqs
	sb.Clear()
	identical = make([][]string, 0)

	// Stores the index of the identical sequence group
	// of a given sequence in the "identical" slice of slice
	seqs := make(map[string]int)
	for _, seq := range oldseqs {
		s := string(seq.sequence)
		// If the group does not exist
		if i, ok := seqs[s]; !ok {
			if err = sb.AddSequence(seq.name, s, seq.comment); err != nil {
				return
			}
			identical = append(identical, []string{seq.name})
			seqs[s] = len(identical) - 1
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

// Returns true if:
//
// - sb and comp have the same number of sequences &&
// - each sequence in sb have a sequence in comp with the same name
//   and the same sequence
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

func (sb *seqbag) IterateChar(it func(name string, sequence []rune) bool) {
	var stop bool = false
	for _, seq := range sb.seqs {
		if stop = it(seq.name, seq.sequence); stop {
			return
		}
	}
}

func (sb *seqbag) IterateAll(it func(name string, sequence []rune, comment string) bool) {
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
func (sb *seqbag) appendToSequence(name string, sequence []rune) error {
	seq, ok := sb.seqmap[name]
	if !ok {
		return fmt.Errorf("Sequence with name %s does not exist in alignment", name)
	}
	seq.sequence = append(seq.sequence, sequence...)
	return nil
}

func (sb *seqbag) AutoAlphabet() {
	isaa := true
	isnt := true

	sb.IterateChar(func(name string, seq []rune) bool {
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
		return false
	})

	if isnt {
		sb.alphabet = NUCLEOTIDS
	} else if isaa {
		sb.alphabet = AMINOACIDS
	} else {
		sb.alphabet = UNKNOWN
	}
}

// Returns the distribution of all characters
func (sb *seqbag) CharStats() map[rune]int64 {
	outmap := make(map[rune]int64)

	for _, seq := range sb.seqs {
		for _, r := range seq.sequence {
			outmap[unicode.ToUpper(r)]++
		}
	}

	return outmap
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
		err = fmt.Errorf("Number of sequences to sample %d is >= sum of the counts %d", nb, total)
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
				err = fmt.Errorf("No sequence named %s is present in the tmp count map", k)
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
	sb.IterateAll(func(name string, sequence []rune, comment string) bool {
		if _, ok := selected[name]; ok {
			sample.AddSequenceChar(name, sequence, comment)
		}
		return false
	})

	return
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
			newseq := []rune(r.ReplaceAllString(string(sb.seqs[seq].sequence), new))
			sb.seqs[seq].sequence = newseq
		}
	} else {
		for seq := 0; seq < sb.NbSequences(); seq++ {
			newseq := strings.Replace(string(sb.seqs[seq].sequence), old, new, -1)
			sb.seqs[seq].sequence = []rune(newseq)
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
		s, _ := sb.seqmap[n]
		sb.seqs[i] = s
	}
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
	var code map[string]rune

	if code, err = geneticCode(geneticcode); err != nil {
		return
	}

	if sb.Alphabet() != NUCLEOTIDS {
		err = errors.New("Wrong alphabet, cannot translate to AA")
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
		err = fmt.Errorf("No ORF has been found on any sequence")
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
		return fmt.Errorf("New name size (%d) does not allow to identify that amount of sequences (%d)",
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
			newlength := int(math.Ceil(math.Log10(float64(*curid + 1))))
			if newlength > length {
				length = newlength
			}
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
