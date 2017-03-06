package align

import (
	"bytes"
	"errors"
	"fmt"
	"github.com/fredericlemoine/goalign/io"
	"math"
	"math/rand"
	"sort"
	"strings"
)

const (
	AMINOACIDS = 0 // Amino acid sequence alphabet
	NUCLEOTIDS = 1 // Nucleotid sequence alphabet
	GAP        = '-'
)

var stdaminoacid = []rune{'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'}
var stdnucleotides = []rune{'A', 'C', 'G', 'T'}

type Alignment interface {
	AddSequence(name string, sequence string, comment string) error
	AddSequenceChar(name string, sequence []rune, comment string) error
	GetSequence(name string) (string, bool)
	GetSequenceChar(ith int) ([]rune, bool)
	GetSequenceName(ith int) (string, bool)
	Iterate(it func(name string, sequence string))
	IterateChar(it func(name string, sequence []rune))
	IterateAll(it func(name string, sequence []rune, comment string))
	NbSequences() int
	Length() int
	ShuffleSequences()
	ShuffleSites(rate float64)
	SimulateRogue(prop float64) ([]string, []string)
	RemoveGaps(cutoff float64)
	Sample(nb int) (Alignment, error)
	BuildBootstrap() Alignment
	Swap(rate float64)
	Recombine(rate float64, lenprop float64)
	Rename(namemap map[string]string)
	TrimNames(size int) (map[string]string, error)
	TrimSequences(trimsize int, fromStart bool) error
	AppendSeqIdentifier(identifier string, right bool)
	AvgAllelesPerSite() float64
	CharStats() map[rune]int64
	Alphabet() int
	Clone() (Alignment, error)
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

func (a *align) IterateAll(it func(name string, sequence []rune, comment string)) {
	for _, seq := range a.seqs {
		it(seq.name, seq.sequence, seq.comment)
	}
}

// Adds a sequence to this alignment
func (a *align) AddSequence(name string, sequence string, comment string) error {
	err := a.AddSequenceChar(name, []rune(sequence), comment)
	return err
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

// If ith >=0 && i < nbSequences() return sequence,true
// Otherwise, return nil, false
func (a *align) GetSequenceChar(ith int) ([]rune, bool) {
	if ith >= 0 && ith < a.NbSequences() {
		return a.seqs[ith].SequenceChar(), true
	}
	return nil, false
}

// If ith >=0 && i < nbSequences() return name,true
// Otherwise, return "", false
func (a *align) GetSequenceName(ith int) (string, bool) {
	if ith >= 0 && ith < a.NbSequences() {
		return a.seqs[ith].Name(), true
	}
	return "", false
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

// Removes positions constituted of [cutoff*100%,100%] Gaps
// Exception fo a cutoff of 0: does not remove positions with 0% gaps
// Cutoff must be between 0 and 1, otherwise set to 0.
// 0 means that positions with > 0 gaps will be removed
// other cutoffs : ]0,1] mean that positions with >= cutoff gaps will be removed
func (a *align) RemoveGaps(cutoff float64) {
	var nbgaps int

	if cutoff < 0 || cutoff > 1 {
		cutoff = 0
	}

	toremove := make([]int, 0, 10)
	for site := 0; site < a.Length(); site++ {
		nbgaps = 0
		for seq := 0; seq < a.NbSequences(); seq++ {
			if a.seqs[seq].sequence[site] == GAP {
				nbgaps++
			}
		}
		if (cutoff > 0.0 && float64(nbgaps) >= cutoff*float64(a.NbSequences())) || (cutoff == 0 && nbgaps > 0) {
			toremove = append(toremove, site)
		}
	}
	/* Now we remove gap positions, starting at the end */
	sort.Ints(toremove)
	for i := (len(toremove) - 1); i >= 0; i-- {
		for seq := 0; seq < a.NbSequences(); seq++ {
			a.seqs[seq].sequence = append(a.seqs[seq].sequence[:toremove[i]], a.seqs[seq].sequence[toremove[i]+1:]...)
		}
	}
	a.length -= len(toremove)
}

// This function renames sequences of the alignment based on the map in argument
// If a name in the map does not exist in the alignment, does nothing
// If a sequence in the alignment does not have a name in the map: does nothing
func (a *align) Rename(namemap map[string]string) {
	for seq := 0; seq < a.NbSequences(); seq++ {
		newname, ok := namemap[a.seqs[seq].name]
		if ok {
			a.seqs[seq].name = newname
		}
		// else {
		// 	io.PrintMessage("Sequence " + a.seqs[seq].name + " not present in the map file")
		// }
	}
}

// Swaps a rate of the sequences together
// takes rate/2 seqs and swap a part of them with the other
// rate/2 seqs at a random position
// if rate < 0 : does nothing
// if rate > 1 : does nothing
func (a *align) Swap(rate float64) {
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
		// We take a random position in the sequences and swap both
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

// Recombines a rate of the sequences to another sequences
// takes rate/2 seqs and copy/paste a portion of them to the other
// rate/2 seqs at a random position
// if rate < 0 : does nothing
// if rate > 1 : does nothing
// prop must be <= 0.5 because it will recombine x% of seqs based on other x% of seqs
func (a *align) Recombine(prop float64, lenprop float64) {
	var seq1, seq2 *seq

	if prop < 0 || prop > 0.5 {
		return
	}
	if lenprop < 0 || lenprop > 1 {
		return
	}

	nb := int(prop * float64(a.NbSequences()))
	lentorecomb := int(lenprop * float64(a.Length()))
	permutation := rand.Perm(a.NbSequences())

	// We take a random position in the sequences between min and max
	for i := 0; i < nb; i++ {
		pos := rand.Intn(a.Length() - lentorecomb + 1)
		seq1 = a.seqs[permutation[i]]
		seq2 = a.seqs[permutation[i+nb]]
		for j := pos; j < pos+lentorecomb; j++ {
			seq1.sequence[j] = seq2.sequence[j]
		}
	}
}

// Simulate rogue taxa in the alignment:
// take the proportion prop of sequences as rogue taxa => R
// For each t in R
//   * We shuffle the alignment sites of t
// Output: List of rogue sequence names, and List of intact sequences
func (a *align) SimulateRogue(prop float64) ([]string, []string) {
	var seq *seq

	if prop < 0 || prop > 1.0 {
		return nil, nil
	}

	nb := int(prop * float64(a.NbSequences()))
	permutation := rand.Perm(a.NbSequences())
	seqlist := make([]string, nb)
	intactlist := make([]string, a.NbSequences()-nb)
	// For each chosen rogue sequence
	for r := 0; r < nb; r++ {
		seq = a.seqs[permutation[r]]
		seqlist[r] = seq.name
		// we Shuffle sequence sites
		for i, _ := range seq.sequence {
			j := rand.Intn(i + 1)
			seq.sequence[i], seq.sequence[j] = seq.sequence[j], seq.sequence[i]
		}
	}
	for nr := nb; nr < a.NbSequences(); nr++ {
		seq = a.seqs[permutation[nr]]
		intactlist[nr-nb] = seq.name
	}
	return seqlist, intactlist
}

func (a *align) TrimNames(size int) (map[string]string, error) {
	shortmap := make(map[string][]string)
	finalshort := make(map[string]string)
	if math.Pow10(size-2) < float64(a.NbSequences()) {
		return nil, errors.New("New name size (" + fmt.Sprintf("%d", size-2) + ") does not allow to identify that amount of sequences (" + fmt.Sprintf("%d", a.NbSequences()) + ")")
	}

	for _, seq := range a.seqs {
		n := seq.Name()
		n = strings.Replace(n, ":", "", -1)
		n = strings.Replace(n, "_", "", -1)
		// Possible to have 100 Identical shortnames
		if len(n) >= size-2 {
			n = n[0 : size-2]
		} else {
			target := len(n)
			for m := 0; m < (size - 2 - target); m++ {
				n = n + "x"
			}
		}
		if _, ok := shortmap[n]; !ok {
			shortmap[n] = make([]string, 0, 100)
		}
		shortmap[n] = append(shortmap[n], seq.Name())
	}

	for short, list := range shortmap {
		if len(list) > 100 {
			return nil, errors.New("More than 100 identical short names (" + short + "), cannot shorten the names")
		} else if len(list) > 1 {
			i := 1
			for _, longname := range list {
				finalshort[longname] = short + fmt.Sprintf("%02d", i)
				i++
			}
		} else {
			for _, longname := range list {
				finalshort[longname] = short + "00"
			}
		}
	}

	for long, short := range finalshort {
		if seq, ok := a.seqmap[long]; !ok {
			return nil, errors.New("The sequence with name " + long + " does not exist")
		} else {
			seq.name = short
			delete(a.seqmap, long)
			a.seqmap[short] = seq
		}
	}

	return finalshort, nil
}

// Trims alignment sequences.
// If fromStart, then trims from the start, else trims from the end
// If trimsize >= sequence or trimsize < 0 lengths, then throw an error
func (a *align) TrimSequences(trimsize int, fromStart bool) error {
	if trimsize < 0 {
		return errors.New("Trim size must not be < 0")
	}
	if trimsize >= a.Length() {
		return errors.New("Trim size must be < alignment length (" + fmt.Sprintf("%d", a.Length()) + ")")
	}
	for _, seq := range a.seqs {
		if fromStart {
			seq.sequence = seq.sequence[trimsize:len(seq.sequence)]
		} else {
			seq.sequence = seq.sequence[0 : len(seq.sequence)-trimsize]
		}
	}
	a.length = a.length - trimsize
	return nil
}

// Append a string to all sequence names of the alignment
// If right is true, then append it to the right of each names,
// otherwise, appends it to the left
func (a *align) AppendSeqIdentifier(identifier string, right bool) {
	if len(identifier) != 0 {
		for _, seq := range a.seqs {
			if right {
				seq.name = seq.name + identifier
			} else {
				seq.name = identifier + seq.name
			}
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

	for _, nt := range strings.ToUpper(seq) {
		couldbent := false
		couldbeaa := false
		switch nt {
		case 'A', 'C', 'B', 'R', 'G', '?', '-', 'D', 'K', 'S', 'H', 'M', 'N', 'V', 'X', 'T', 'W', 'Y':
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
			io.ExitWithMessage(errors.New("Unknown character state in alignment : " + string(nt)))
		}
	}

	if isnt {
		return NUCLEOTIDS
	} else {
		return AMINOACIDS
	}
}

// Returns the distribution of all characters
func (a *align) CharStats() map[rune]int64 {
	outmap := make(map[rune]int64)

	for _, seq := range a.seqs {
		for _, r := range seq.sequence {
			outmap[r]++
		}
	}

	return outmap
}

func (a *align) Alphabet() int {
	return a.alphabet
}

func RandomAlignment(alphabet, length, nbseq int) (Alignment, error) {
	al := NewAlign(alphabet)
	for i := 0; i < nbseq; i++ {
		name := fmt.Sprintf("Seq%04d", i)
		if seq, err := RandomSequence(alphabet, length); err != nil {
			return nil, err
		} else {
			al.AddSequenceChar(name, seq, "")
		}
	}
	return al, nil
}

func (a *align) Clone() (Alignment, error) {
	c := NewAlign(a.Alphabet())
	var err error
	a.IterateAll(func(name string, sequence []rune, comment string) {
		newseq := make([]rune, 0, len(sequence))
		newseq = append(newseq, sequence...)
		err = c.AddSequenceChar(name, newseq, comment)
		if err != nil {
			return
		}
	})
	return c, err
}

func (a *align) AvgAllelesPerSite() float64 {
	nballeles := 0
	nbsites := 0
	for site := 0; site < a.Length(); site++ {
		alleles := make(map[rune]bool)
		onlygap := true
		for seq := 0; seq < a.NbSequences(); seq++ {
			s := a.seqs[seq].sequence[site]
			if s != GAP {
				alleles[s] = true
				onlygap = false
			}
		}
		if !onlygap {
			nbsites++
		}
		nballeles += len(alleles)
	}
	return float64(nballeles) / float64(nbsites)
}
