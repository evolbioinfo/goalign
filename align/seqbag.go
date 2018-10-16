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

	"github.com/fredericlemoine/goalign/io"
)

type SeqBag interface {
	AddSequence(name string, sequence string, comment string) error
	AddSequenceChar(name string, sequence []rune, comment string) error
	AppendSeqIdentifier(identifier string, right bool)
	Alphabet() int
	AlphabetStr() string
	AlphabetCharacters() []rune
	AutoAlphabet() // detects and sets alphabet automatically for all the sequences
	CharStats() map[rune]int64
	CleanNames()                            // Clean sequence names (newick special char)
	Clear()                                 // Removes all sequences
	Deduplicate() error                     // Remove duplicate sequences
	GetSequence(name string) (string, bool) // Get a sequence by names
	GetSequenceById(ith int) (string, bool)
	GetSequenceChar(name string) ([]rune, bool)
	GetSequenceCharById(ith int) ([]rune, bool)
	GetSequenceNameById(ith int) (string, bool)
	SetSequenceChar(ithAlign, ithSite int, char rune) error
	Sequence(ith int) (Sequence, bool)
	SequenceByName(name string) (Sequence, bool)
	Identical(SeqBag) bool
	Iterate(it func(name string, sequence string))
	IterateChar(it func(name string, sequence []rune))
	IterateAll(it func(name string, sequence []rune, comment string))
	LongestORF() (orf Sequence, err error)
	NbSequences() int
	Phase(orf Sequence) (seqs SeqBag, aaseqs SeqBag, positions []int, err error)
	Rename(namemap map[string]string)
	RenameRegexp(regex, replace string, namemap map[string]string) error
	ShuffleSequences()               // Shuffle sequence order
	Translate(phase int) (err error) // Translates nt sequence in aa
	TrimNames(namemap map[string]string, size int) error
	TrimNamesAuto(namemap map[string]string, curid *int) error
	Sort() // Sorts the sequences by name
	Unalign() SeqBag
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

// Removes spaces and tabs at beginning and end of sequence names
// and replaces newick special characters \s\t()[];,.: by "-"
func (sb *seqbag) CleanNames() {
	firstlast := regexp.MustCompile("(^[\\s\\t]+|[\\s\\t]+$)")
	inside := regexp.MustCompile("[\\s\\t,\\[\\]\\(\\),;\\.:]+")

	for _, seq := range sb.seqs {
		seq.name = firstlast.ReplaceAllString(seq.name, "")
		seq.name = inside.ReplaceAllString(seq.name, "-")
	}
}

// Removes all the sequences from the seqbag
func (sb *seqbag) Clear() {
	sb.seqmap = make(map[string]*seq)
	sb.seqs = make([]*seq, 0, 100)
}

// This function removes sequences that are duplicates of other
// It keeps one copy of each sequence, with the name of the first
// found.
//
// It modifies input alignment.
func (sb *seqbag) Deduplicate() (err error) {
	oldseqs := sb.seqs
	sb.Clear()

	seqs := make(map[string]bool)
	for _, seq := range oldseqs {
		s := string(seq.sequence)
		_, ok := seqs[s]
		if !ok {
			if err = sb.AddSequence(seq.name, s, seq.comment); err != nil {
				return
			}
			seqs[s] = true
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

	if isnt {
		return NUCLEOTIDS
	} else {
		return AMINOACIDS
	}
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
- if the alphabet is not NUCLEOTIDES: returns an error

The seqbag is cleared and old sequences are replaced with aminoacid sequences

It only translates using the standard code so far.
*/
func (sb *seqbag) Translate(phase int) (err error) {
	var oldseqs []*seq
	var buffer bytes.Buffer

	if sb.Alphabet() != NUCLEOTIDS {
		err = errors.New("Wrong alphabet, cannot translate to AA")
		return
	}

	oldseqs = sb.seqs
	sb.Clear()

	for _, seq := range oldseqs {
		buffer.Reset()

		if len(seq.sequence) < 3+phase {
			err = fmt.Errorf("Cannot translate a sequence with length < 3+phase (%s)", seq.name)
			return
		}
		for i := phase; i < len(seq.sequence)-2; i += 3 {
			codon := strings.Replace(strings.ToUpper(string(seq.sequence[i:i+3])), "U", "T", -1)
			aa, found := standardcode[codon]
			if !found {
				aa = 'X'
			}
			buffer.WriteRune(aa)
		}
		if err2 := sb.AddSequence(seq.name, buffer.String(), seq.comment); err != nil {
			err = err2
			return
		}
	}

	return
}

func (sb *seqbag) LongestORF() (orf Sequence, err error) {
	var beststart, bestend int
	var start, end int
	var bestseq Sequence

	// Search for the longest orf in all sequences
	for _, seq := range sb.seqs {
		start, end = seq.LongestORF()
		if end-start > bestend-beststart {
			beststart, bestend = start, end
			bestseq = seq
		}
		if start == -1 {
			err = fmt.Errorf("No ORF is found on the sequence %s", seq.name)
			return
		}
	}
	// log.Print("Longest ORF found in sequence ", bestseq.Name())
	// log.Print(string(bestseq.SequenceChar()[beststart:bestend]))
	orf = NewSequence("longestorf_sw", bestseq.SequenceChar()[beststart:bestend], "")
	return
}

// align all sequences to the given ORF and trims sequences to the start
// position
// If orf is nil, searches for the longest ORF in all sequences
//
// To do so, Phase() will:
//
// 1. Translate the given ORF in aminoacids;
// 2. For each sequence of the dataset: translate it in the 3 phases (forward strand only),
//    align it with the translated orf, and take the phase giving the best alignment; If no phase
//    gives a good alignment (>80% orf length, and starting at first position of the ORF), then
//    the sequence is discarded;
// 3. For each sequence, take the Start corresponding to the Start of the ORF, and remove
//    nucleotides before;
// 4. Return the trimmed nucleotidic sequences (phased), the corresponding amino-acid sequences (phasedaa)
//    and the positions of starts in the nucleotidic sequences.
//
// It does not modify the input object
func (sb *seqbag) Phase(orf Sequence) (phased SeqBag, phasedaa SeqBag, positions []int, err error) {
	var aligner PairwiseAligner
	var orfaa, seqaa, bestseq Sequence
	var phase, beststart, beststartaa int
	var bestscore float64

	if sb.Alphabet() != NUCLEOTIDS {
		err = fmt.Errorf("Wrong alphabet for phase : %s", sb.AlphabetStr())
		return
	}
	positions = make([]int, 0, sb.NbSequences())
	phased = NewSeqBag(sb.Alphabet())
	phasedaa = NewSeqBag(AMINOACIDS)

	if orf == nil {
		orf, err = sb.LongestORF()
	}

	// We translate the longest ORF in AA
	if orfaa, err = orf.Translate(0); err != nil {
		return
	}

	// Now we align all sequences against this longest orf aa sequence with Modified Smith/Waterman
	for _, seq := range sb.seqs {
		bestscore = 0.0
		beststart = 0
		beststartaa = 0
		bestseq = nil
		// We translate the sequence in the 3 phases to search for the best
		// alignment
		for phase = 0; phase < 3; phase++ {
			if seqaa, err = seq.Translate(phase); err != nil {
				return
			}
			aligner = NewPwAligner(orfaa, seqaa, ALIGN_ALGO_ATG)
			aligner.SetGapScore(-4.0)
			aligner.SetMatchScore(2.0)
			aligner.SetMismatchScore(-4.0)

			aligner.Alignment()
			orfstart, seqstart := aligner.AlignStarts()
			_, seqend := aligner.AlignEnds()

			if orfstart == 0 && float64(seqend-seqstart)/float64(orfaa.Length()) >= 0.8 {
				if aligner.MaxScore() > bestscore {
					bestscore = aligner.MaxScore()
					// Alignment start in nucleotidic sequence
					beststart = phase + (seqstart * 3)
					// Alignment start in proteic sequence
					beststartaa = seqstart
					bestseq = seqaa
				}
			}
		}

		if bestscore > 0.0 {
			positions = append(positions, beststart)
			phased.AddSequence(seq.name, string(seq.sequence[beststart:]), seq.comment)
			phasedaa.AddSequence(seq.name, string(bestseq.SequenceChar()[beststartaa:]), seq.comment)
		} else {
			log.Print(seq.name, "\n", "Cannot find a good ORF alignment for sequence", "\n")
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
