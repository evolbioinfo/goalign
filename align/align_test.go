package align

import (
	"fmt"
	"math"
	"strings"
	"testing"
)

func TestRandomAlignment(t *testing.T) {
	length := 3000
	nbseqs := 500
	a, err := RandomAlignment(AMINOACIDS, length, nbseqs)
	if err != nil {
		t.Error(err)
	}

	if a.Length() != length {
		t.Error(fmt.Sprintf("Length should be %d and is %d", length, a.Length()))
	}
	if a.NbSequences() != nbseqs {
		t.Error(fmt.Sprintf("Nb sequences should be %d and is %d", nbseqs, a.NbSequences()))
	}
}

func TestAppendIdentifier(t *testing.T) {
	a, err := RandomAlignment(AMINOACIDS, 300, 50)
	if err != nil {
		t.Error(err)

	}
	a.AppendSeqIdentifier("IDENT", false)

	a.IterateChar(func(name string, sequence []rune) {
		if !strings.HasPrefix(name, "IDENT") {
			t.Error("Sequence name does not start with expected id: IDENT")
		}
	})

	a.AppendSeqIdentifier("IDENT", true)
	a.IterateChar(func(name string, sequence []rune) {
		if !strings.HasSuffix(name, "IDENT") {
			t.Error("Sequence name does not end with expected id: IDENT")
		}
	})
}

func TestCleanNames(t *testing.T) {
	a, err := RandomAlignment(AMINOACIDS, 300, 50)
	if err != nil {
		t.Error(err)

	}
	a2, _ := a.Clone()
	a.AppendSeqIdentifier("\t \t", false)
	a.AppendSeqIdentifier("\t \t", true)

	a.CleanNames()
	i := 0
	a.IterateChar(func(name string, sequence []rune) {
		expected, found := a2.GetSequenceNameById(i)
		if !found {
			t.Error("Unknown sequence name after clean")
		}
		if name != expected {
			t.Error("Unexpected sequence name after clean")
		}
		i++
	})
}

func TestRemoveOneGapSite(t *testing.T) {
	a, err := RandomAlignment(AMINOACIDS, 300, 300)
	if err != nil {
		t.Error(err)

	}

	/* We add 1 gap per site */
	pos := 0
	a.IterateChar(func(name string, sequence []rune) {
		sequence[pos] = GAP
		pos++
	})

	a.RemoveGapSites(0.0)

	if a.Length() != 0 {
		t.Error("We should have removed all positions")
	}
	a.IterateChar(func(name string, sequence []rune) {
		if len(sequence) != 0 {
			t.Error(fmt.Sprintf("Sequence length after removing gaps should be 0 and is : %d", len(sequence)))
		}
	})
}

func TestRemoveAllGapSites(t *testing.T) {
	a, err := RandomAlignment(AMINOACIDS, 300, 300)
	if err != nil {
		t.Error(err)

	}

	backupseq := make([]rune, 0, 300)
	seq0, found := a.GetSequenceCharById(0)
	if !found {
		t.Error("Problem finding first sequence")
	}

	/* We add all gaps on 1 site */
	/* And one gap at all sites */
	pos1 := 20
	pos2 := 0
	a.IterateChar(func(name string, sequence []rune) {
		sequence[pos1] = GAP
		sequence[pos2] = GAP
		pos2++
	})
	backupseq = append(backupseq, seq0...)
	/* Remove position 20 */
	backupseq = append(backupseq[:20], backupseq[21:]...)

	a.RemoveGapSites(1.0)

	if a.Length() != 299 {
		t.Error("We should have removed only one position")
	}

	a.IterateChar(func(name string, sequence []rune) {
		if len(sequence) != 299 {
			t.Error(fmt.Sprintf("Sequence length after removing gaps should be 299 and is : %d", len(sequence)))
		}
	})

	newseq, found2 := a.GetSequenceCharById(0)
	if !found2 {
		t.Error("Problem finding first seqence")
	}

	for i, c := range newseq {
		if c != backupseq[i] {
			t.Error(fmt.Sprintf("Char at position %d should be %c and is %c", i, backupseq[i], c))
		}
	}
}

func TestRemoveOneGapSequence(t *testing.T) {
	a, err := RandomAlignment(AMINOACIDS, 300, 300)
	if err != nil {
		t.Error(err)

	}

	/* We add 1 gap per site on all sequences*/
	pos := 0
	a.IterateChar(func(name string, sequence []rune) {
		sequence[pos] = GAP
		pos++
	})

	a.RemoveGapSeqs(0.0)

	if a.NbSequences() != 0 {
		t.Error("We should have removed all sequences")
	}
}

func TestRemoveOneGapSequence2(t *testing.T) {
	a, err := RandomAlignment(AMINOACIDS, 300, 300)
	if err != nil {
		t.Error(err)

	}

	/* We add 1 gap per site on half of the sequences*/
	pos := 0
	a.IterateChar(func(name string, sequence []rune) {
		if pos%2 == 0 {
			sequence[pos] = GAP
		}
		pos++
	})

	a.RemoveGapSeqs(0.0)

	if a.NbSequences() != 150 {
		t.Error("We should have removed half of sequences")
	}
}

func TestRemoveAllGapSequences(t *testing.T) {
	a, err := RandomAlignment(AMINOACIDS, 300, 300)
	if err != nil {
		t.Error(err)

	}
	seq0, found := a.GetSequenceCharById(0)
	if !found {
		t.Error("Problem finding first sequence")
	}

	for i := 0; i < a.Length(); i++ {
		seq0[i] = GAP
	}

	a.RemoveGapSeqs(1.0)

	if a.NbSequences() != 299 {
		t.Error("We should have removed only one sequence")
	}
}

func TestRemoveHalfGapSequences(t *testing.T) {
	a, err := RandomAlignment(AMINOACIDS, 300, 300)
	if err != nil {
		t.Error(err)

	}
	seq0, found := a.GetSequenceCharById(0)
	if !found {
		t.Error("Problem finding first sequence")
	}

	for i := 0; i < a.Length(); i++ {
		if i%2 == 0 {
			seq0[i] = GAP
		}
	}

	a.RemoveGapSeqs(0.5)

	if a.NbSequences() != 299 {
		t.Error("We should have removed only one sequence")
	}
}

func TestClone(t *testing.T) {
	a, err := RandomAlignment(AMINOACIDS, 300, 300)
	if err != nil {
		t.Error(err)

	}

	/* We add 1 gap per site */
	pos := 0
	a.IterateChar(func(name string, sequence []rune) {
		sequence[pos] = GAP
		pos++
	})

	a2, err2 := a.Clone()
	if err2 != nil {
		t.Error(err2)
	}

	a.RemoveGapSites(0.0)

	a2.IterateChar(func(name string, sequence []rune) {
		if len(sequence) != 300 {
			t.Error(fmt.Sprintf("Clone lenght should be 300 and is : %d", len(sequence)))
		}
	})
}

func TestClone2(t *testing.T) {
	a, err := RandomAlignment(AMINOACIDS, 300, 300)
	if err != nil {
		t.Error(err)

	}

	a2, err2 := a.Clone()
	if err2 != nil {
		t.Error(err2)
	}

	i := 0
	a2.IterateChar(func(name string, sequence []rune) {
		s2, ok := a.GetSequenceCharById(i)
		n2, ok2 := a.GetSequenceNameById(i)

		if !ok || !ok2 {
			t.Error(fmt.Sprintf("Sequence not found in clone alignment: %s", name))
		}

		if len(sequence) != len(s2) {
			t.Error(fmt.Sprintf("Clone length is different from original length : %d != %d", len(sequence), len(s2)))
		}
		if name != n2 {
			t.Error(fmt.Sprintf("Clone and original sequences at position %d have different names : %s != %s", i, name, n2))
		}
		for j, c := range sequence {
			if c != s2[j] {
				t.Error(fmt.Sprintf("Clone sequence is different from original at position %d : %c != %c", j, c, s2[j]))
			}
		}
		i++
	})
}

func TestAvgAlleles(t *testing.T) {
	a, err := RandomAlignment(AMINOACIDS, 300, 300)
	if err != nil {
		t.Error(err)

	}

	a.IterateChar(func(name string, sequence []rune) {
		for j, _ := range sequence {
			sequence[j] = 'A'
		}
	})

	if a.AvgAllelesPerSite() != 1 {
		t.Error("There should be 1 allele per site in this alignment")
	}
}

func TestAvgAlleles2(t *testing.T) {
	a, err := RandomAlignment(AMINOACIDS, 300, 300)
	if err != nil {
		t.Error(err)

	}

	i := 0
	a.IterateChar(func(name string, sequence []rune) {
		for j, _ := range sequence {
			/* One only gap sequence */
			if i == 10 {
				sequence[j] = GAP
			} else if i <= 75 {
				sequence[j] = 'A'
			} else if i <= 150 {
				sequence[j] = 'C'
			} else if i <= 225 {
				sequence[j] = 'G'
			} else {
				sequence[j] = 'T'
			}
		}
		// Add a gap at a whole position
		sequence[50] = GAP
		i++
	})
	fmt.Println(a.AvgAllelesPerSite())
	if a.AvgAllelesPerSite() != 4 {
		t.Error("There should be 4 allele per site in this alignment")
	}
}

func TestRename(t *testing.T) {
	length := 3000
	nbseqs := 500
	a, err := RandomAlignment(AMINOACIDS, length, nbseqs)

	if err != nil {
		t.Error(err)
	}

	// Map File
	namemap1 := make(map[string]string)
	namemap2 := make(map[string]string)
	for i := 0; i < 500; i++ {
		namemap1[fmt.Sprintf("Seq%04d", i)] = fmt.Sprintf("New%04d", i)
		namemap2[fmt.Sprintf("New%04d", i)] = fmt.Sprintf("Seq%04d", i)
	}

	a.Rename(namemap1)
	i := 0
	a.IterateChar(func(name string, sequence []rune) {
		expected := fmt.Sprintf("New%04d", i)
		if name != expected {
			t.Error(fmt.Sprintf("Sequence name should be %s and is %s", expected, name))
		}
		i++
	})

	a.Rename(namemap2)
	i = 0
	a.IterateChar(func(name string, sequence []rune) {
		expected := fmt.Sprintf("Seq%04d", i)
		if name != expected {
			t.Error(fmt.Sprintf("Sequence name should be %s and is %s", expected, name))
		}
		i++
	})
}

func TestRogue(t *testing.T) {
	length := 3000
	nbseqs := 500
	nrogue := 0.5
	a, err := RandomAlignment(AMINOACIDS, length, nbseqs)

	if err != nil {
		t.Error(err)
	}

	a2, err2 := a.Clone()
	if err2 != nil {
		t.Error(err2)
	}

	rogues, intacts := a2.SimulateRogue(nrogue, 1.0)

	if (len(rogues) + len(intacts)) != a.NbSequences() {
		t.Error("Number of intact + rogue sequences is not the same than the total number of sequences")
	}

	for _, s := range rogues {
		seq2, ok2 := a2.GetSequence(s)
		seq, ok := a.GetSequence(s)
		if !ok || !ok2 {
			t.Error("Rogue name does not exist in alignment...")
		}
		if seq == seq2 {
			t.Error("Rogue sequence is the same after simulation (should be shuffled)...")
		}
	}
	for _, s := range intacts {
		seq2, ok2 := a2.GetSequence(s)
		seq, ok := a.GetSequence(s)
		if !ok || !ok2 {
			t.Error("Intact name does not exist in alignment...")
		}
		if seq != seq2 {
			t.Error("Intact sequences should be he same before and after simulation...")
		}
	}
}

func TestRogue2(t *testing.T) {
	length := 3000
	nbseqs := 500
	nrogue := 0.0
	a, err := RandomAlignment(AMINOACIDS, length, nbseqs)

	if err != nil {
		t.Error(err)
	}

	a2, err2 := a.Clone()
	if err2 != nil {
		t.Error(err2)
	}

	rogues, intacts := a2.SimulateRogue(nrogue, 1.0)

	if (len(rogues) + len(intacts)) != a.NbSequences() {
		t.Error("Number of intact + rogue sequences is not the same than the total number of sequences")
	}

	if len(rogues) > 0 {
		t.Error("There should be no rogue taxa in output")
	}

	for _, s := range intacts {
		seq2, ok2 := a2.GetSequence(s)
		seq, ok := a.GetSequence(s)
		if !ok || !ok2 {
			t.Error("Intact name does not exist in alignment...")
		}
		if seq != seq2 {
			t.Error("Intact sequences should be he same before and after simulation...")
		}
	}
}

func TestRogue3(t *testing.T) {
	length := 3000
	nbseqs := 500
	nrogue := 1.0
	a, err := RandomAlignment(AMINOACIDS, length, nbseqs)

	if err != nil {
		t.Error(err)
	}

	a2, err2 := a.Clone()
	if err2 != nil {
		t.Error(err2)
	}

	rogues, intacts := a2.SimulateRogue(nrogue, 1.0)

	if (len(rogues) + len(intacts)) != a.NbSequences() {
		t.Error(fmt.Sprintf("Number of intact (%d) + rogue (%d) sequences is not the same than the total number of sequences (%d)", len(intacts), len(rogues), a.NbSequences()))
	}

	if len(rogues) < a.NbSequences() {
		t.Error("All sequences should be rogues")
	}

	if len(intacts) > 0 {
		t.Error("There should be no intact sequences")
	}
	for _, s := range rogues {
		seq2, ok2 := a2.GetSequence(s)
		seq, ok := a.GetSequence(s)
		if !ok || !ok2 {
			t.Error("Rogue name does not exist in alignment...")
		}
		if seq == seq2 {
			t.Error("Rogue sequence is the same after simulation (should be shuffled)...")
		}
	}
}

func TestEntropy(t *testing.T) {
	length := 3
	nbseqs := 5
	a, err := RandomAlignment(AMINOACIDS, length, nbseqs)

	alldifferent := []rune{'A', 'R', 'N', 'D', 'C'}
	// First site: only 'R' => Entropy 0.0
	// Second site: Only different aminoacids => Entropy 1.0
	for i := 0; i < 5; i++ {
		a.SetSequenceChar(i, 0, 'R')
		a.SetSequenceChar(i, 1, alldifferent[i])
	}
	e, err := a.Entropy(0, false)
	if err != nil {
		t.Error(err)
	}
	if e != 0.0 {
		t.Error(fmt.Sprintf("Entropy should be 0.0 and is %f", e))
	}

	e, err = a.Entropy(1, false)
	if err != nil {
		t.Error(err)
	}
	expected := -1.0 * float64(nbseqs) * 1.0 / float64(nbseqs) * math.Log(1.0/float64(nbseqs))
	if int(e*10000000) != int(expected*10000000) {
		t.Error(fmt.Sprintf("Entropy should be %.7f and is %.7f", expected, e))
	}
}

func TestSubAlign(t *testing.T) {
	length := 200
	nbseqs := 1001
	a, err := RandomAlignment(AMINOACIDS, length, nbseqs)

	/* We put only A from index 9 to index 99 */
	for i := 0; i < 1001; i++ {
		for aa := 9; aa < 99; aa++ {
			a.SetSequenceChar(i, aa, 'A')
		}
	}

	subalign, err := a.SubAlign(9, 90)
	if err != nil {
		t.Error(err)
	}
	subalign.IterateChar(func(name string, sequence []rune) {
		if len(sequence) != 90 {
			t.Error(fmt.Sprintf("Length of subsequence must be %d and is %d", 90, len(sequence)))
		}
		for i, a := range sequence {
			if a != 'A' {
				t.Error(fmt.Sprintf("Character at position %d must be %c and is %c", i, 'A', a))
			}
		}
	})
}

func TestConcat(t *testing.T) {
	var err error
	var a Alignment
	var a2 Alignment
	var acopy Alignment
	length := 200
	nbseqs := 1001
	a, err = RandomAlignment(AMINOACIDS, length, nbseqs)
	if err != nil {
		t.Error(err)
	}
	a2, err = RandomAlignment(AMINOACIDS, length*2, nbseqs)
	if err != nil {
		t.Error(err)
	}

	acopy, err = a.Clone()
	if err != nil {
		t.Error(err)
	}

	err = a.Concat(a2)
	if err != nil {
		t.Error(err)
	}

	if a.Length() != length*3 {
		t.Error(fmt.Sprintf("Concatenated alignment should have a length of %d and not %d", length*3, a.Length()))
	}

	a.IterateChar(func(name string, sequence []rune) {
		s, _ := acopy.GetSequence(name)
		s2, _ := a2.GetSequence(name)
		if string(sequence) != s+s2 {
			t.Error(fmt.Sprintf("Concatenated sequence is not correct"))
		}
	})
}

func TestDedup(t *testing.T) {
	var err error
	var a Alignment

	a = NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "ACGT", "")
	a.AddSequence("B", "ACGG", "")
	a.AddSequence("C", "ACGT", "")
	a.AddSequence("C", "ACGT", "")
	a.AddSequence("D", "ACGT", "")

	if a.NbSequences() != 5 {
		t.Error("There should be 5 sequences before deduplicaion")
	}

	if err = a.Deduplicate(); err != nil {
		t.Error(err)
	}

	if a.NbSequences() != 2 {
		t.Error("There should be 2 sequences in the deduplicated alignment")
	}

	if err = a.Deduplicate(); err != nil {
		t.Error(err)
	}
	if a.NbSequences() != 2 {
		t.Error("There should be 2 sequences in the de-deduplicated alignment")
	}
}

func TestCodonAlign(t *testing.T) {
	a := NewAlign(UNKNOWN)
	a.AddSequence("Seq0000", "D*-AVGQNLK", "")
	a.AddSequence("Seq0001", "IE-FKF-LLM", "")
	a.AddSequence("Seq0002", "ERTSSYFLNT", "")
	a.AutoAlphabet()

	n := NewSeqBag(UNKNOWN)
	n.AddSequence("Seq0000", "GATTAAGCCGTAGGCCAGAATCTGAAG", "")
	n.AddSequence("Seq0001", "ATCGAATTTAAGTTTCTTCTAATG", "")
	n.AddSequence("Seq0002", "GAGAGGACTAGTTCATACTTTTTAAACACT", "")
	n.AutoAlphabet()

	exp := NewSeqBag(UNKNOWN)
	exp.AddSequence("Seq0000", "GATTAA---GCCGTAGGCCAGAATCTGAAG", "")
	exp.AddSequence("Seq0001", "ATCGAA---TTTAAGTTT---CTTCTAATG", "")
	exp.AddSequence("Seq0002", "GAGAGGACTAGTTCATACTTTTTAAACACT", "")
	exp.AutoAlphabet()

	res, err := a.CodonAlign(n)
	if err != nil {
		t.Error(err)
	}

	if res.Alphabet() != NUCLEOTIDS {
		t.Error(fmt.Errorf("Alphabet of result codon align must be nucleotides"))
	}
	if a.Alphabet() != AMINOACIDS {
		t.Error(fmt.Errorf("Alphabet of aa alignment must be aminoacids"))
	}
	if n.Alphabet() != NUCLEOTIDS {
		t.Error(fmt.Errorf("Alphabet of nt seqs to align must be nucleotides"))
	}
	if !exp.Identical(res) {
		t.Error(fmt.Errorf("Expected alignment is different from codon aligned alignemnts"))
	}
}

func TestIdenticalAligns(t *testing.T) {
	var err error
	var a, a2, a3 Alignment

	if a, err = RandomAlignment(AMINOACIDS, 1001, 145); err != nil {
		t.Error(err)
	}

	if a2, err = a.Clone(); err != nil {
		t.Error(err)
	}

	if a3, err = a.Clone(); err != nil {
		t.Error(err)
	}

	a3.ShuffleSequences()

	if !a.Identical(a) {
		t.Error(fmt.Errorf("Same alignment object must be Identical"))
	}

	if !a.Identical(a2) {
		t.Error(fmt.Errorf("Cloned alignment must be Identical"))
	}

	if !a.Identical(a3) {
		t.Error(fmt.Errorf("Shuffled alignment must be Identical"))
	}
}
