package align

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

type phasedsequence struct {
	removed  bool
	position int
	// phased nt sequence
	ntseq *seq
	// phased aa sequence
	aaseq *seq

	ali *align

	// positions of frameshifts in the sequence
	// compared to reference sequence/orf phase
	frameshifts []int
	// positions of stop codons in the sequence
	// compared to reference sequence/orf phase.
	// Only stop codons before end of sequence.
	stops []int
}

func (p *phaser) Phase(orfs SeqBag, seqs SeqBag) (phased <-chan phasedsequence, err error) {
	if p.translate {
		phased, err = p.phase(orfs, seqs)
	} else {
		phased, err = p.phasent(orfs, seqs)
	}
	return
}

// orfs: Reference sequences/ORFs to phase sequences with
// seqs: Sequences to phase
func (p *phaser) phase(orfs SeqBag, seqs SeqBag) (phased <-chan phasedsequence, err error) {
	phased = make(chan phasedsequence, 50)

	var orf Sequence
	var alphabet int
	var lock, lock2 sync.Mutex

	// Channels for concurrency
	seqchan := make(chan *seq)
	if seqs.Alphabet() != NUCLEOTIDS {
		err = fmt.Errorf("Wrong alphabet for phase : %s", sb.AlphabetStr())
		return
	}

	if orfs == nil {
		if orf, err = sb.LongestORF(p.reverse); err != nil {
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
	go func() {
		for _, seq := range seqs.Sequences() {
			seqchan <- seq
		}
		close(seqchan)
	}()

	// All threads consuming sequences
	var wg sync.WaitGroup
	for cpu := 0; cpu < cpus; cpu++ {
		wg.Add(1)
		go func(cpu int) {
			var bestscore float64
			var bestratematches, bestlen float64
			var beststart, bestend int
			var beststartaa, bestendaa int
			var bestseq, bestseqaa Sequence
			var phase int
			var phases int // Number of phases 3 or 6
			var seqaa Sequence
			var tmpseq, revcomp Sequence
			var aligner PairwiseAligner
			var p *phasedsequence
			var al *align

			defer wg.Done()

			for seq := range seqchan {
				if err != nil {
					return
				}
				bestscore = .0
				beststart = 0
				bestend = 0
				beststartaa = 0
				bestendaa = 0
				bestseq = nil
				bestratematches = .0
				bestlen = .0

				// We translate the sequence in the 3 phases to search for the best
				// alignment
				phases = 3
				tmpseq = seq
				revcomp = seq
				if reverse {
					phases = 6
					revcomp = seq.Clone()
					revcomp.Reverse()
					revcomp.Complement()
				}

				// We search for the best score among all references
				// and all phases
				for _, orfaa := range orfsaa.Sequences() {
					for phase = 0; phase < phases; phase++ {
						if phase < 3 {
							tmpseq = seq
						} else {
							tmpseq = revcomp
						}
						if seqaa, err = tmpseq.Translate(phase % 3); err != nil {
							return
						}
						aligner = NewPwAligner(orfaa, seqaa, ALIGN_ALGO_ATG)
						aligner.SetGapOpenScore(-10.0)
						aligner.SetGapExtendScore(-.5)

						if al, err = aligner.Alignment(); err != nil {
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
							bestseq2ali = a.Seq2Ali()
							bestratematches = float64(aligner.NbMatches()) / float64(aligner.Length())
							bestlen = float64(aligner.Length()) / float64(orfaa.Length())
							bestend = bestseq.Length()
							bestendaa = bestseqaa.Length()
							if cutend {
								bestend = (phase % 3) + ((seqend + 1) * 3)
								bestendaa = seqend + 1
							}
						}
					}
				}

				p = phasedsequence{
					removed:  false,
					position: beststart,
					ntseq: NewSequence(bestseq.Name(),
						bestseq.SequenceChar()[beststart:bestend],
						bestseq.Comment()),
					aaseq: NewSequence(bestseqaa.Name(),
						string(bestseqaa.SequenceChar()[beststartaa:bestendaa]),
						bestseqaa.Comment()),
					ali: bestali,
				} // We set a threshold for matches over the alignment length...
				if (matchcutoff > .0 && bestratematches <= matchcutoff) ||
					(lencutoff > .0 && bestlen <= lencutoff) {
					p.removed = true
				}
				phased <- p
			}
		}(cpu)
	}
	wg.Wait()
	return
}

// orfs: Reference sequences/ORFs to phase sequences with
// seqs: Sequences to phase
func (p *phaser) phasent(orfs SeqBag, seqs SeqBag) (phased <-chan phasedsequence, err error) {
	phased = make(chan phasedsequence, 50)

	var orf Sequence
	var alphabet int
	var lock, lock2 sync.Mutex

	// Channels for concurrency
	seqchan := make(chan *seq)
	if seqs.Alphabet() != NUCLEOTIDS {
		err = fmt.Errorf("Wrong alphabet for phasent : %s", sb.AlphabetStr())
		return
	}

	if orfs == nil {
		if orf, err = sb.LongestORF(p.reverse); err != nil {
			return
		}
		orfs = NewSeqBag(UNKNOWN)
		orfs.AddSequenceChar(orf.Name(), orf.SequenceChar(), orf.Comment())
		orfs.AutoAlphabet()
	}

	// We translate the longest ORF in AA if it is nucleotides
	alphabet = orfs.Alphabet()
	if alphabet != NUCLEOTIDS && alphabet != BOTH {
		err = fmt.Errorf("Wrong orf alphabet")
		return
	}

	// Now we align all sequences against this longest orf aa sequence with Modified Smith/Waterman
	// We use n threads
	// Fill the sequence channel
	go func() {
		for _, seq := range seqs.Sequences() {
			seqchan <- seq
		}
		close(seqchan)
	}()

	// All threads consuming sequences
	var wg sync.WaitGroup
	for cpu := 0; cpu < cpus; cpu++ {
		wg.Add(1)
		go func(cpu int) {
			var bestscore float64
			var bestratematches, bestlen float64
			var beststart, bestend int
			var bestseq Sequence
			var phase int
			var phases int // Number of phases 1 or 2 (+/+-)
			var tmpseq, revcomp Sequence
			var aligner PairwiseAligner
			var p *phasedsequence
			var bestali *align

			defer wg.Done()

			for seq := range seqchan {
				if err != nil {
					return
				}
				bestscore = .0
				beststart = 0
				bestend = 0
				bestseq = nil
				bestratematches = .0
				bestlen = .0

				// We translate the sequence in the 3 phases to search for the best
				// alignment
				phases = 1
				tmpseq = seq
				revcomp = seq
				if reverse {
					phases = 2
					revcomp = seq.Clone()
					revcomp.Reverse()
					revcomp.Complement()
				}
				// We search for the best score among all references
				// and all phases
				for _, orfaa := range orfsaa.Sequences() {
					for phase = 0; phase < phases; phase++ {
						if phase < 2 {
							tmpseq = seq
						} else {
							tmpseq = revcomp
						}
						aligner = NewPwAligner(orf, tmpseq, ALIGN_ALGO_ATG)
						aligner.SetGapOpenScore(-10.0)
						aligner.SetGapExtendScore(-.5)

						if al, err = aligner.Alignment(); err != nil {
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
							bestratematches = float64(aligner.NbMatches()) / float64(aligner.Length())
							bestlen = float64(aligner.Length()) / float64(orfaa.Length())
							bestend = bestseq.Length()
							if cutend {
								bestend = seqend + 1
							}
						}
					}
				}

				ph = &phasedsequence{
					removed:  false,
					position: beststart,
					ntseq: NewSequence(bestseq.Name(),
						bestseq.SequenceChar()[beststart:bestend],
						bestseq.Comment()),
					aaseq: NewSequence(bestseqaa.Name(),
						string(bestseqaa.SequenceChar()[beststartaa:bestendaa]),
						bestseqaa.Comment()),
					ali: al,
				} // We set a threshold for matches over the alignment length...
				if err = ph.aaseq.Translate(0); err != nil {
					return
				}
				if (matchcutoff > .0 && bestratematches <= matchcutoff) ||
					(lencutoff > .0 && bestlen <= lencutoff) {
					ph.removed = true
				} else {
					ph.computestops()
				}
				phased <- p
			}
		}(cpu)
	}
	wg.Wait()
	return
}
