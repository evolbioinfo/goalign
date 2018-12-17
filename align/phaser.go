package align

import (
	"fmt"
	"sync"
)

// * If SetTranslate(true):
//
// align all sequences to the given ORF and trims sequences to the start
// position
// If orf is nil, searches for the longest ORF (in 3 or 6 phases depending on reverse arg) in all sequences
//
// To do so, Phase() will:
//
// 1. Translate the given ORF in aminoacids;
// 2. For each sequence of the dataset: translate it in the 3 phases (forward) if reverse is false or 6
//    phases (forward and reverse) if reverse is true, align it with the translated orf, and take the phase
//    giving the best alignment; If no phase gives a good alignment (>lencutoff * orf length, >matchcutoff
//    matches over the align length and starting at first position of the ORF), then the sequence is discarded;
// 3. For each sequence, take the Start corresponding to the Start of the ORF, and remove
//    nucleotides before;
// 4. Return the trimmed nucleotidic sequences (phased), the corresponding amino-acid sequences (phasedaa)
//    the positions of starts in the nucleotidic sequences, and the removed sequence names.
//
// If cutend is true, then also remove the end of sequences that do not align with orf
//
// It does not modify the input object
//
//
// * If SetTranslate(false):
//
// align all sequences to the given ORF and trims sequences to the start
// position, it does not take into account protein information
//
// If orf is nil, searches for the longest ORF (in forward only or both strands depending on reverse arg) in all sequences
//
// To do so:
//
// 1. If alignment is bad (>lencutoff * orf length, >matchcutoff matches over the align length and starting at first position of the ORF), then the sequence is discarded;
// 3. For each sequence, take the Start corresponding to the Start of the ORF, and remove nucleotides before;
// 4. Return the trimmed nucleotidic sequences (phased), the positions of starts in the nucleotidic sequences, and the removed sequence names.
// If cutend is true, then also remove the end of sequences that do not align with orf
// It does not modify the input object
type Phaser interface {
	Phase(orfs, seqs SeqBag) (chan PhasedSequence, error)
	SetLenCutoff(cutoff float64)
	SetMatchCutoff(cutoff float64)
	SetReverse(reverse bool)
	SetCutEnd(cutend bool)
	SetCpus(cpus int)
	SetTranslate(translate bool)
}

type phaser struct {
	// Alignment length cutoff for sequences to orfs alignment
	lencutoff float64
	// #Matches cutoff for sequences to orfs alignment
	matchcutoff float64
	// If sequences must be analyzed also in reverse strand
	reverse bool
	// Cut the end of sequences
	cutend bool
	// Number of CPUs for computation
	cpus int
	// Translate orf (1st phase) & input sequences in 3 or 6 phases
	// to search for the best phase. If false, just align sequences
	// as is
	translate bool
}

type PhasedSequence struct {
	Err      error
	Removed  bool
	Position int
	// phased nt sequence
	NtSeq Sequence
	// phased aa sequence
	AaSeq Sequence
	// Aligned sequences
	// 1st: best found orf
	// 2nd: sequence
	Ali Alignment
}

func NewPhaser() Phaser {
	return &phaser{
		lencutoff:   .8,
		matchcutoff: .5,
		reverse:     false,
		cutend:      false,
		cpus:        1,
		translate:   true,
	}
}

func (p *phaser) SetLenCutoff(cutoff float64) {
	p.lencutoff = cutoff
}
func (p *phaser) SetMatchCutoff(cutoff float64) {
	p.matchcutoff = cutoff
}
func (p *phaser) SetReverse(reverse bool) {
	p.reverse = reverse
}
func (p *phaser) SetCutEnd(cutend bool) {
	p.cutend = cutend
}
func (p *phaser) SetCpus(cpus int) {
	p.cpus = cpus
}
func (p *phaser) SetTranslate(translate bool) {
	p.translate = translate
}

// orfs: Reference sequences/ORFs to phase sequences with
// seqs: Sequences to phase
func (p *phaser) Phase(orfs, seqs SeqBag) (phased chan PhasedSequence, err error) {
	phased = make(chan PhasedSequence, 50)

	var orf Sequence
	var alphabet int
	var orfsaa SeqBag

	// Channels for concurrency
	var seqchan <-chan Sequence
	if seqs.Alphabet() != NUCLEOTIDS {
		err = fmt.Errorf("Wrong alphabet for phase : %s", seqs.AlphabetStr())
		return
	}

	// If no orf given, then we find the longest among the sequences
	if orfs == nil {
		if orf, err = seqs.LongestORF(p.reverse); err != nil {
			return
		}
		orfs = NewSeqBag(UNKNOWN)
		orfs.AddSequenceChar(orf.Name(), orf.SequenceChar(), orf.Comment())
		orfs.AutoAlphabet()
	}

	// We translate the longest ORF in AA if it is nucleotides
	alphabet = orfs.Alphabet()
	if orfsaa, err = orfs.CloneSeqBag(); err != nil {
		return
	}
	if alphabet == NUCLEOTIDS {
		if err = orfsaa.Translate(0); err != nil {
			return
		}
	}

	// Now we align all sequences against this longest orf aa sequence with Modified Smith/Waterman
	// We use n threads
	// Fill the sequence channel
	seqchan = seqs.SequencesChan()

	// All threads consuming sequences
	var wg sync.WaitGroup
	var goerr error
	for cpu := 0; cpu < p.cpus; cpu++ {
		wg.Add(1)
		go func() {
			var inerr error
			defer wg.Done()
			var ph PhasedSequence

			for seq := range seqchan {
				if p.translate {
					ph, inerr = p.alignAgainstRefsAA(seq, orfsaa.Sequences())
				} else {
					ph, inerr = p.alignAgainstRefsNT(seq, orfs.Sequences())
				}
				if inerr != nil {
					goerr = inerr
				}
				if goerr != nil {
					return
				}
				phased <- ph
			}
		}()
	}

	go func() {
		wg.Wait()
		close(phased)
		// In case an error occured
		// we must finish to read the seqchan
		for _ = range seqchan {
		}
	}()
	return
}

func (p *phaser) alignAgainstRefsAA(seq Sequence, orfsaa []Sequence) (ph PhasedSequence, err error) {
	var bestscore float64 = .0
	var bestratematches, bestlen float64 = .0, .0
	var beststart, bestend int = 0, 0
	var beststartaa, bestendaa int = 0, 0
	var bestseq, bestseqaa Sequence = nil, nil
	var bestali Alignment = nil
	var phases int = 3 // Number of phases 3 or 6

	var phase int
	var seqaa Sequence
	var tmpseq, revcomp Sequence
	var aligner PairwiseAligner
	var al Alignment

	// We translate the sequence in the 3 phases to search for the best
	// alignment
	tmpseq = seq
	revcomp = seq
	if p.reverse {
		phases = 6
		revcomp = seq.Clone()
		revcomp.Reverse()
		revcomp.Complement()
	}

	// We search for the best score among all references
	// and all phases
	for _, orfaa := range orfsaa {
		for phase = 0; phase < phases; phase++ {
			if phase < 3 {
				tmpseq = seq
			} else {
				tmpseq = revcomp
			}
			if seqaa, err = tmpseq.Translate(phase % 3); err != nil {
				ph = PhasedSequence{Err: fmt.Errorf("Error while translating %s : %v", seq.Name(), err)}
				return
			}
			aligner = NewPwAligner(orfaa, seqaa, ALIGN_ALGO_ATG)
			aligner.SetGapOpenScore(-10.0)
			aligner.SetGapExtendScore(-.5)

			if al, err = aligner.Alignment(); err != nil {
				ph = PhasedSequence{Err: fmt.Errorf("Error while aligning %s with %s : %v", orfaa.Name(), seqaa.Name(), err)}
				return
			}
			_, seqstart := aligner.AlignStarts()
			_, seqend := aligner.AlignEnds()

			if aligner.MaxScore() > bestscore {
				bestscore = aligner.MaxScore()
				// Alignment start in nucleotidic sequence
				beststart = (phase % 3) + (seqstart * 3)
				// Alignment start in proteic sequence
				beststartaa = seqstart
				bestseqaa = seqaa
				bestseq = tmpseq
				bestali = al
				bestratematches = float64(aligner.NbMatches()) / float64(aligner.Length())
				bestlen = float64(aligner.Length()) / float64(orfaa.Length())
				bestend = bestseq.Length()
				bestendaa = bestseqaa.Length()
				if p.cutend {
					bestend = (phase % 3) + ((seqend + 1) * 3)
					bestendaa = seqend + 1
				}
			}
		}
	}

	ph = PhasedSequence{
		Err:      nil,
		Removed:  false,
		Position: beststart,
		NtSeq: NewSequence(bestseq.Name(),
			bestseq.SequenceChar()[beststart:bestend],
			bestseq.Comment()),
		AaSeq: NewSequence(bestseqaa.Name(),
			bestseqaa.SequenceChar()[beststartaa:bestendaa],
			bestseqaa.Comment()),
		Ali: bestali,
	} // We set a threshold for matches over the alignment length...
	if (p.matchcutoff > .0 && bestratematches <= p.matchcutoff) ||
		(p.lencutoff > .0 && bestlen <= p.lencutoff) {
		ph.Removed = true
	}
	return
}

// orfs: Reference sequences/ORFs to phase sequences with
// seqs: Sequences to phase
func (p *phaser) alignAgainstRefsNT(seq Sequence, orfs []Sequence) (ph PhasedSequence, err error) {
	var bestscore float64 = .0
	var bestratematches, bestlen float64 = .0, .0
	var beststart, bestend int = 0, 0
	var bestseq Sequence = nil
	var phases int = 1     // Number of phases 1 or 2 (+/+-)
	var nbgapstart int = 0 // Number of gaps at the start of compared seq
	var phase int
	var tmpseq, revcomp Sequence
	var aligner PairwiseAligner
	var al, bestali Alignment

	// We translate the sequence in the 3 phases to search for the best
	// alignment
	phases = 1
	tmpseq = seq
	revcomp = seq
	if p.reverse {
		phases = 2
		revcomp = seq.Clone()
		revcomp.Reverse()
		revcomp.Complement()
	}

	// We search for the best score among all references
	// and all phases
	for _, orf := range orfs {
		for phase = 0; phase < phases; phase++ {
			if phase < 1 {
				tmpseq = seq
			} else {
				tmpseq = revcomp
			}
			aligner = NewPwAligner(orf, tmpseq, ALIGN_ALGO_ATG)
			aligner.SetGapOpenScore(-10.0)
			aligner.SetGapExtendScore(-.5)

			if al, err = aligner.Alignment(); err != nil {
				ph = PhasedSequence{Err: fmt.Errorf("Error while aligning %s with %s : %v", orf.Name(), tmpseq.Name(), err)}
				return
			}
			_, seqstart := aligner.AlignStarts()
			_, seqend := aligner.AlignEnds()

			if aligner.MaxScore() > bestscore {
				bestscore = aligner.MaxScore()
				// Alignment start in nucleotidic sequence
				beststart = seqstart
				bestseq = tmpseq
				bestali = al
				// Number of gaps at the start of sequence
				// Allows to translate in the right phase
				nbgapstart = 0
				for i := 0; aligner.Seq2Ali()[i] == '-'; i++ {
					nbgapstart++
				}
				bestratematches = float64(aligner.NbMatches()) / float64(aligner.Length())
				bestlen = float64(aligner.Length()) / float64(orf.Length())
				bestend = bestseq.Length()
				if p.cutend {
					bestend = seqend + 1
				}
			}
		}
	}

	ph = PhasedSequence{
		Err:      nil,
		Removed:  false,
		Position: beststart,
		NtSeq: NewSequence(bestseq.Name(),
			bestseq.SequenceChar()[beststart:bestend],
			bestseq.Comment()),
		AaSeq: NewSequence(bestseq.Name(),
			bestseq.SequenceChar()[beststart:bestend],
			bestseq.Comment()),
		Ali: bestali,
	}

	// We translate in the right phase (taking into account initial gaps)
	// Ex:
	// NNN NNN NNN => phase 0
	// -NN NNN NNN => phase 2
	// --N NNN NNN => phase 1
	// --- NNN NNN => phase 0
	ph.AaSeq, err = ph.AaSeq.Translate((3 - nbgapstart%3) % 3)

	// We set a threshold for matches over the alignment length...
	if (p.matchcutoff > .0 && bestratematches <= p.matchcutoff) ||
		(p.lencutoff > .0 && bestlen <= p.lencutoff) {
		ph.Removed = true
	}
	return
}
