package align

import (
	"fmt"
	"math"
	"sort"
	"strings"
	"testing"
)

func TestRanTranslatedomAlignment(t *testing.T) {
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

	a.CleanNames(nil)
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

func TestReplace(t *testing.T) {
	length := 3000
	nbseqs := 500
	a, err := RandomAlignment(NUCLEOTIDS, length, nbseqs)

	if err != nil {
		t.Error(err)
	}

	acount := 0
	a.IterateChar(func(name string, sequence []rune) {
		for _, c := range sequence {
			if c == 'A' {
				acount++
			}
		}
	})

	gapcount := 0
	a.Replace("A", "-", false)
	a.IterateChar(func(name string, sequence []rune) {
		for _, c := range sequence {
			if c == 'A' {
				t.Error(fmt.Sprintf("There should not remains A after replace"))
			}
			if c == '-' {
				gapcount++
			}
		}
	})
	if gapcount != acount {
		t.Error(fmt.Sprintf("Each A should have been replaced by a -"))
	}
}

func TestReplaceRegEx(t *testing.T) {
	length := 3000
	nbseqs := 500
	a, err := RandomAlignment(NUCLEOTIDS, length, nbseqs)

	if err != nil {
		t.Error(err)
	}

	regcount := 0
	a.IterateChar(func(name string, sequence []rune) {
		prev := '0'
		for _, c := range sequence {
			if (c == 'A' && prev == 'A') ||
				(c == 'C' && prev == 'A') ||
				(c == 'G' && prev == 'A') ||
				(c == 'T' && prev == 'A') {
				regcount += 2
				prev = '-'
			} else {
				prev = c
			}
		}
	})

	gapcount := 0
	a.Replace("A.", "--", true)
	a.IterateChar(func(name string, sequence []rune) {
		prev := '0'
		for _, c := range sequence {
			if (c == 'A' && prev == 'A') ||
				(c == 'C' && prev == 'A') ||
				(c == 'G' && prev == 'A') ||
				(c == 'T' && prev == 'A') {
				t.Error(fmt.Sprintf("There should not remains %c%c after replace", prev, c))
			}
			if c == '-' {
				gapcount++
			}
			prev = c
		}
	})

	if gapcount != regcount {
		t.Error(fmt.Sprintf("Each A. should have been replaced by -- (%d != %d)", gapcount, regcount))
	}
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

func TestUnalign(t *testing.T) {

	in := NewAlign(UNKNOWN)
	in.AddSequence("Seq0000", "GATTAA---GCCGTAGGCCAGAATCTGAAG", "")
	in.AddSequence("Seq0001", "ATCGAA---TTTAAGTTT---CTTCTAATG", "")
	in.AddSequence("Seq0002", "GAGAGGACTAGTTCATACTTTTTAAACACT", "")
	in.AutoAlphabet()

	exp := NewSeqBag(UNKNOWN)
	exp.AddSequence("Seq0000", "GATTAAGCCGTAGGCCAGAATCTGAAG", "")
	exp.AddSequence("Seq0001", "ATCGAATTTAAGTTTCTTCTAATG", "")
	exp.AddSequence("Seq0002", "GAGAGGACTAGTTCATACTTTTTAAACACT", "")
	exp.AutoAlphabet()

	res := in.Unalign()

	if !exp.Identical(res) {
		t.Error(fmt.Errorf("Expected sequences are different from unaligned alignemnt"))
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

func TestPhase(t *testing.T) {

	in := NewSeqBag(UNKNOWN)
	in.AddSequence("Seq0000", "GCTGCCGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGATAGC", "")
	in.AddSequence("Seq0001", "GCTGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA", "")
	in.AutoAlphabet()

	exp := NewSeqBag(UNKNOWN)
	exp.AddSequence("Seq0000", "ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGATAGC", "")
	exp.AddSequence("Seq0001", "ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA", "")
	exp.AutoAlphabet()

	expaa := NewSeqBag(UNKNOWN)
	expaa.AddSequence("Seq0000", "MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVVV***", "")
	expaa.AddSequence("Seq0001", "MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVV**", "")
	expaa.AutoAlphabet()

	phaser := NewPhaser()
	phaser.SetLenCutoff(-1.0)
	phaser.SetMatchCutoff(0.5)
	phaser.SetReverse(true)
	phaser.SetCutEnd(false)
	phaser.SetCpus(1)

	seqs := NewSeqBag(UNKNOWN)
	seqsaa := NewSeqBag(UNKNOWN)

	phased, err := phaser.Phase(nil, in)
	if err != nil {
		t.Error(err)
	}

	for ph := range phased {
		if ph.Err != nil {
			t.Error(ph.Err)
		}
		seqs.AddSequence(ph.NtSeq.Name(), ph.NtSeq.Sequence(), ph.NtSeq.Comment())
		seqsaa.AddSequence(ph.AaSeq.Name(), ph.AaSeq.Sequence(), ph.AaSeq.Comment())
	}
	seqs.AutoAlphabet()
	seqsaa.AutoAlphabet()
	if !exp.Identical(seqs) {
		t.Error(fmt.Errorf("Expected sequences are different from phased sequences \n%s \n%s", exp.String(), seqs.String()))
	}
	if !expaa.Identical(seqsaa) {
		t.Error(fmt.Errorf("Expected sequences are different from phased sequences"))
	}
}

func TestPhaseReverse(t *testing.T) {

	in := NewSeqBag(UNKNOWN)
	in.AddSequence("Seq0000", "GCTATCATTACACTACGACAACTATGATAATGTAATAGTGATGCCACCCTCCGCCACCCGTTGTGGTAGTCTCTTCGCTACTCGATGAGGAAGACTGTTGCGGTGGGGGAGGGCAACAGAAAAAGTCATCCATGTTATTCTTTTTCCTTCTCCGTCGGCGACGCAGTAGGAGAAGCAATAACGCTGCGGCAGC", "")
	in.AddSequence("Seq0001", "GCTGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA", "")
	in.AutoAlphabet()

	exp := NewSeqBag(UNKNOWN)
	exp.AddSequence("Seq0000", "ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGATAGC", "")
	exp.AddSequence("Seq0001", "ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA", "")
	exp.AutoAlphabet()

	expaa := NewSeqBag(UNKNOWN)
	expaa.AddSequence("Seq0000", "MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVVV***", "")
	expaa.AddSequence("Seq0001", "MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVV**", "")
	expaa.AutoAlphabet()

	phaser := NewPhaser()
	phaser.SetLenCutoff(-1.0)
	phaser.SetMatchCutoff(0.5)
	phaser.SetReverse(true)
	phaser.SetCutEnd(false)
	phaser.SetCpus(1)

	seqs := NewSeqBag(UNKNOWN)
	seqsaa := NewSeqBag(UNKNOWN)

	phased, err := phaser.Phase(nil, in)
	if err != nil {
		t.Error(err)
	}

	for ph := range phased {
		if ph.Err != nil {
			t.Error(ph.Err)
		}
		seqs.AddSequence(ph.NtSeq.Name(), ph.NtSeq.Sequence(), ph.NtSeq.Comment())
		seqsaa.AddSequence(ph.AaSeq.Name(), ph.AaSeq.Sequence(), ph.AaSeq.Comment())
	}
	seqs.AutoAlphabet()
	seqsaa.AutoAlphabet()

	if !exp.Identical(seqs) {
		t.Error(fmt.Errorf("Expected sequences are different from phased sequences (with reverse) \n%s \n%s", exp.String(), seqs.String()))
	}
	if !expaa.Identical(seqsaa) {
		t.Error(fmt.Errorf("Expected sequences are different from phased sequences (with reverse)"))
	}
}

func TestPhaseReverseCutEnds(t *testing.T) {

	in := NewSeqBag(UNKNOWN)
	in.AddSequence("Seq0000", "GCTATCATTACACTACGACAACTATGATAATGTAATAGTGATGCCACCCTCCGCCACCCGTTGTGGTAGTCTCTTCGCTACTCGATGAGGAAGACTGTTGCGGTGGGGGAGGGCAACAGAAAAAGTCATCCATGTTATTCTTTTTCCTTCTCCGTCGGCGACGCAGTAGGAGAAGCAATAACGCTGCGGCAGC", "")
	in.AddSequence("Seq0001", "GCTGCAGCGTTATTGCTTCTCCTACTGCGTCGCCGACGGAGAAGGAAAAAGAATAACATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA", "")
	in.AutoAlphabet()

	exp := NewSeqBag(UNKNOWN)
	exp.AddSequence("Seq0000", "ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAA", "")
	exp.AddSequence("Seq0001", "ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA", "")
	exp.AutoAlphabet()

	expaa := NewSeqBag(UNKNOWN)
	expaa.AddSequence("Seq0000", "MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVVV*", "")
	expaa.AddSequence("Seq0001", "MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVV**", "")
	expaa.AutoAlphabet()

	phaser := NewPhaser()
	phaser.SetLenCutoff(-1.0)
	phaser.SetMatchCutoff(0.5)
	phaser.SetReverse(true)
	phaser.SetCutEnd(true)
	phaser.SetCpus(1)

	seqs := NewSeqBag(UNKNOWN)
	seqsaa := NewSeqBag(UNKNOWN)

	phased, err := phaser.Phase(nil, in)
	if err != nil {
		t.Error(err)
	}

	for ph := range phased {
		if ph.Err != nil {
			t.Error(ph.Err)
		}
		seqs.AddSequence(ph.NtSeq.Name(), ph.NtSeq.Sequence(), ph.NtSeq.Comment())
		seqsaa.AddSequence(ph.AaSeq.Name(), ph.AaSeq.Sequence(), ph.AaSeq.Comment())
	}
	seqs.AutoAlphabet()
	seqsaa.AutoAlphabet()

	if !exp.Identical(seqs) {
		t.Error(fmt.Errorf("Expected sequences are different from phased sequences (with reverse) \n%s \n%s", exp.String(), seqs.String()))
	}
	if !expaa.Identical(seqsaa) {
		t.Error(fmt.Errorf("Expected sequences are different from phased sequences (with reverse) \n%s \n%s", expaa.String(), seqsaa.String()))
	}
}

func TestTranslate(t *testing.T) {

	in := NewSeqBag(UNKNOWN)
	in.AddSequence("Seq0000", "ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGATAGC", "")
	in.AddSequence("Seq0001", "ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA", "")
	in.AutoAlphabet()

	in2, err2 := in.CloneSeqBag()
	if err2 != nil {
		t.Error(err2)
	}

	expaa := NewSeqBag(UNKNOWN)
	expaa.AddSequence("Seq0000", "MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVVV***", "")
	expaa.AddSequence("Seq0001", "MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVV**", "")
	expaa.AutoAlphabet()

	exp3phases := NewSeqBag(UNKNOWN)
	exp3phases.AddSequence("Seq0000_0", "MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVVV***", "")
	exp3phases.AddSequence("Seq0000_1", "WMTFSVALPHRNSLPHRVAKRLPQRVAEGGITITLS*LS*CNDS", "")
	exp3phases.AddSequence("Seq0000_2", "G*LFLLPSPTATVFLIE*RRDYHNGWRRVASLLHYHSCRSVMI", "")
	exp3phases.AddSequence("Seq0001_0", "MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIIVVV**", "")
	exp3phases.AddSequence("Seq0001_1", "WMTFSVALPHRNSLPHRVAKRLPQRVAEGGITITLS*LSYN", "")
	exp3phases.AddSequence("Seq0001_2", "G*LFLLPSPTATVFLIE*RRDYHNGWRRVASLLHYHSCRIM", "")
	exp3phases.AutoAlphabet()

	if err := in.Translate(0, 0); err != nil {
		t.Error(err)
	} else {
		if !expaa.Identical(in) {
			t.Error(fmt.Errorf("Expected sequences are different from phased sequences"))
		}
	}

	if err := in2.Translate(-1, 0); err != nil {
		t.Error(err)
	} else {
		if !exp3phases.Identical(in2) {
			t.Error(fmt.Errorf("Expected sequences are different from 3 phase translated sequences"))
		}
	}
}

func TestTranslateMito(t *testing.T) {

	in := NewSeqBag(UNKNOWN)
	in.AddSequence("Seq0000", "ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTAGTGTAATGATAGC", "")
	in.AddSequence("Seq0001", "ATGGATGACTTTTTCTGTTGCCCTCCCCCACCGCAACAGTCTTCCTCATCGAGTAGCGAAGAGACTACCACAACGGGTGGCGGAGGGTGGCATCACTATTACATTATCATAGTTGTCGTATAATGA", "")
	in.AutoAlphabet()

	in2, err2 := in.CloneSeqBag()
	if err2 != nil {
		t.Error(err2)
	}

	expaa := NewSeqBag(UNKNOWN)
	expaa.AddSequence("Seq0000", "MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIMVVVV*W*", "")
	expaa.AddSequence("Seq0001", "MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIMVVV*W", "")
	expaa.AutoAlphabet()

	exp3phases := NewSeqBag(UNKNOWN)
	exp3phases.AddSequence("Seq0000_0", "MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIMVVVV*W*", "")
	exp3phases.AddSequence("Seq0000_1", "WMTFSVALPHRNSLPHRVAK*LPQRVAEGGITITLS*LS*CNDS", "")
	exp3phases.AddSequence("Seq0000_2", "GWLFLLPSPTATVFLIE*R*DYHNGWR*VASLLHYHSCRSVMM", "")
	exp3phases.AddSequence("Seq0001_0", "MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIMVVV*W", "")
	exp3phases.AddSequence("Seq0001_1", "WMTFSVALPHRNSLPHRVAK*LPQRVAEGGITITLS*LSYN", "")
	exp3phases.AddSequence("Seq0001_2", "GWLFLLPSPTATVFLIE*R*DYHNGWR*VASLLHYHSCRMM", "")
	exp3phases.AutoAlphabet()

	if err := in.Translate(0, 1); err != nil {
		t.Error(err)
	} else {
		if !expaa.Identical(in) {
			t.Error(fmt.Errorf("Expected sequences are different from phased sequences"))
		}
	}

	if err := in2.Translate(-1, 1); err != nil {
		t.Error(err)
	} else {
		if !exp3phases.Identical(in2) {
			t.Error(fmt.Errorf("Expected sequences are different from 3 phase translated sequences"))
		}
	}
}

func TestMaskNt(t *testing.T) {
	in := NewAlign(UNKNOWN)
	in.AddSequence("Seq0000", "GATTAATTTGCCGTAGGCCA", "")
	in.AddSequence("Seq0001", "GAATCTGAAGATCGAACACT", "")
	in.AddSequence("Seq0002", "TTAAGTTTTCACTTCTAATG", "")
	in.AutoAlphabet()

	exp := NewSeqBag(UNKNOWN)
	exp.AddSequence("Seq0000", "GATTAATTNNCCGTAGGCCA", "")
	exp.AddSequence("Seq0001", "GAATCTGANNATCGAACACT", "")
	exp.AddSequence("Seq0002", "TTAAGTTTNNACTTCTAATG", "")
	exp.AutoAlphabet()

	if err := in.Mask(8, 2); err != nil {
		t.Error(err)
	} else {
		if !exp.Identical(in) {
			t.Error(fmt.Errorf("Expected sequences are different from masked sequences"))
		}
	}
}

func TestMaskProt(t *testing.T) {
	in := NewAlign(UNKNOWN)
	in.AddSequence("Seq0000", "PHGVHCVSSYRFEKCPNFFC", "")
	in.AddSequence("Seq0001", "EACKWDNTCPMKIETHQHQK", "")
	in.AddSequence("Seq0002", "GDMMEDSGSIAIDGIGHHKN", "")
	in.AutoAlphabet()

	exp := NewSeqBag(UNKNOWN)
	exp.AddSequence("Seq0000", "PHGVHCVSSYRFXXXXXXXX", "")
	exp.AddSequence("Seq0001", "EACKWDNTCPMKXXXXXXXX", "")
	exp.AddSequence("Seq0002", "GDMMEDSGSIAIXXXXXXXX", "")
	exp.AutoAlphabet()

	if err := in.Mask(12, 2000); err != nil {
		t.Error(err)
	} else {
		if !exp.Identical(in) {
			t.Error(fmt.Errorf("Expected sequences are different from masked sequences"))
		}
	}
}

func TestDiff(t *testing.T) {
	in := NewAlign(UNKNOWN)
	in.AddSequence("Seq0000", "GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTTTTTAAACACTTTTACATCGA", "")
	in.AddSequence("Seq0001", "TGTCGGACCTAAGTATTGAGTACAACGGTGTATTCCAGCGGTGGAGAGGTCTATTTTTCCGGTTGAAGGACTCTAGAGCTGTAAAGGGTATGGCCATGTG", "")
	in.AddSequence("Seq0002", "CTAAGCGCGGGCGGATTGCTGTTGGAGCAAGGTTAAATACTCGGCAATGCCCCATGATCCCCCAAGGACAATAAGAGCGAAGTTAGAACAAATGAACCCC", "")
	in.AddSequence("Seq0003", "GAGTGGAGGCTTTATGGCACAAGGTATTAGAGACTGAGGGGCACCCCGGCATGGTAAGCAGGAGCCATCGCGAAGGCTTCAGGTATCTTCCTGTGTTACC", "")
	in.AddSequence("Seq0004", "CATAGCCCCTGATGCCCTGACCCGTGTCGCGGCAACGTCTACATTTCACGATAAATACTCCGCTGCTAGTCGGCTCTAGATGCTTTTCTTCCAGATCTGG", "")
	in.AddSequence("Seq0005", "AGTTTGACTATGAGCGCCGGCTTAGTGCTGACAGTGATGCTCCGTTGTAAGGGTCCTGATGTTCTTGTGCTCGCGCATATTAGAGCTGAGTTTCCCAAAG", "")
	in.AddSequence("Seq0006", "TCGCCACGGTGTGGAATGTACGTTATGGCAGTAATCAGCGGCTTTCACCGACATGCCCCCTCCGTGGCTCCTTGCGACCATCGGCGGACCTGCGGTGTCG", "")
	in.AddSequence("Seq0007", "CTGGTAATACCTGCGCTATTTCGTCAGTTCGTGTACGGGTAACGATAGCGGTTAATGCTTATTCCGATCAGCTCACACCCATGAAGGTGGCTCTGGAGCC", "")
	in.AddSequence("Seq0008", "TCGTTAACCCACTCTAACCACCTCCTGTAGCGACATCGGGTGCTCGGCTTGGATACCTTCGTCATATTGGACCCCAGGTCTCAACCTCGTGAGCTCTCTG", "")
	in.AddSequence("Seq0009", "ACCTACGGCTCTAGACAGCTGAAGTCCGGTTCCGAGCACTGTACGGAAACTTGAAAAGGCTCGACGGAGGCTTGTTCCGCAGAGTGGGACTATAACATAC", "")
	in.AutoAlphabet()

	exp := NewAlign(UNKNOWN)
	exp.AddSequence("Seq0000", "GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTTTTTAAACACTTTTACATCGA", "")
	exp.AddSequence("Seq0001", "TG.CGGACCTAA...TTGAGT.CAAC.GT.TATTCCAG.GG.GGAGAGGTCTA.T.TTCC.GTT.A.GG.C.CT.G.GC.G.A..GGGTA.GGC...GTG", "")
	exp.AddSequence("Seq0002", "CTAAGCGCG.G..G.TTG.T.TTGGA.C.AGGTT..ATAC.CGGCAA.G.C.CATG.TCCCCC.A.GAC.A.AAGAG.GAAG.T.GA..AAA.GA.C.CC", "")
	exp.AddSequence("Seq0003", "..G.GGAGGCTTTAT...ACA.GGTATT...GACTGAGGGGC.CCCCGG..TGGTA.GCA.GAGCC.TCGCGAAGGCT.CAGGT.T.TTCC.GTGT.ACC", "")
	exp.AddSequence("Seq0004", "C..AGCCCCTGATGCCCTG.CCCGTGTCGCGG.A.CGT..AC.TT.CACG.TAAA..C.CCGCT.CTAGTCGG.TCTAGA.GCTTTTCT.CCAGATCT.G", "")
	exp.AddSequence("Seq0005", "AG..TGAC.ATGAGC.C.GGCTTAG..CT..CA.TGATGC.CCGT.G.AAGGG..CTGAT.TTCTTGTGCTCG.GC.TA..AG.GCTGAG...C.CAAAG", "")
	exp.AddSequence("Seq0006", "TCGCC.CGGTGT.G.ATGT.CGT.A..GCAG.AATCAG.GGCTTTCACCG..A.GCCCCCTCCGT.G..CC..GCG..CA.CGGCGG..C.GCGGTGTCG", "")
	exp.AddSequence("Seq0007", "CTGGT.A.AC.T.CGCTATTTCG..A.TTCG.GT.CGGG.AACGA.AGCGGT.AA.GC.TATTCC..TC..C...C..CCA.G..GGTGGC.CTGGAGCC", "")
	exp.AddSequence("Seq0008", "TCG.T.ACCCA.TCTAA...CCTC...T..CGAC.T.GGG.GCTCGGC.TGGA.ACCT.C.TC.TATTGGACC.CAGG.C.CA.CCTCG.GAGCTC..TG", "")
	exp.AddSequence("Seq0009", "ACC..CGGCT.TAG.CAG.T...GTCCGGTTC...G....G..C.GAAA.TTGAAA.GGCTC..C.GAGGC..GT.C.GCAGAGTGGGAC.A..ACATAC", "")

	in.DiffWithFirst()
	if !exp.Identical(in) {
		t.Error(fmt.Errorf("Expected sequences are different from Diff sequences"))
	}
}

func TestDiff2(t *testing.T) {
	in := NewAlign(UNKNOWN)
	in.AddSequence("Seq0000", "GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTTTTTAAACACTTTTACATCGA", "")
	in.AddSequence("Seq0001", "GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTTTTTAAACACTTTTACATCGA", "")
	in.AddSequence("Seq0002", "GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTTTTTAAACACTTTTACATCGA", "")
	in.AddSequence("Seq0003", "GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTTTTTAAACACTTTTACATCGA", "")
	in.AddSequence("Seq0004", "GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTTTTTAAACACTTTTACATCGA", "")
	in.AddSequence("Seq0005", "GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTTTTTAAACACTTTTACATCGA", "")
	in.AddSequence("Seq0006", "GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTTTTTAAACACTTTTACATCGA", "")
	in.AddSequence("Seq0007", "GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTTTTTAAACACTTTTACATCGA", "")
	in.AddSequence("Seq0008", "GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTTTTTAAACACTTTTACATCGA", "")
	in.AddSequence("Seq0009", "GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTTTTTAAACACTTTTACATCGA", "")
	in.AutoAlphabet()

	exp := NewAlign(UNKNOWN)
	exp.AddSequence("Seq0000", "GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTTTTTAAACACTTTTACATCGA", "")
	exp.AddSequence("Seq0001", "....................................................................................................", "")
	exp.AddSequence("Seq0002", "....................................................................................................", "")
	exp.AddSequence("Seq0003", "....................................................................................................", "")
	exp.AddSequence("Seq0004", "....................................................................................................", "")
	exp.AddSequence("Seq0005", "....................................................................................................", "")
	exp.AddSequence("Seq0006", "....................................................................................................", "")
	exp.AddSequence("Seq0007", "....................................................................................................", "")
	exp.AddSequence("Seq0008", "....................................................................................................", "")
	exp.AddSequence("Seq0009", "....................................................................................................", "")

	in.DiffWithFirst()
	if !exp.Identical(in) {
		t.Error(fmt.Errorf("Expected sequences are different from Diff sequences"))
	}
}

func TestDiffCount(t *testing.T) {
	in := NewAlign(UNKNOWN)
	in.AddSequence("Seq0000", "GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTTTTTAAACACTTTTACATCGA", "")
	in.AddSequence("Seq0001", "TGTCGGACCTAAGTATTGAGTACAACGGTGTATTCCAGCGGTGGAGAGGTCTATTTTTCCGGTTGAAGGACTCTAGAGCTGTAAAGGGTATGGCCATGTG", "")
	in.AddSequence("Seq0002", "GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTTTTTAAACACTTTGACATCGA", "")
	in.AddSequence("Seq0003", "GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTTTTTAAACACTTTTACATCGA", "")
	in.AutoAlphabet()

	exp := []string{"AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"}
	expdiffs := []map[string]int{
		{
			"AC": 5,
			"AG": 12,
			"AT": 5,
			"CA": 5,
			"CG": 5,
			"CT": 6,
			"GA": 2,
			"GC": 2,
			"GT": 8,
			"TA": 7,
			"TC": 7,
			"TG": 10,
		},
		{
			"TG": 1,
		},
		{},
	}
	res, resdiffs := in.CountDifferences()

	sort.Strings(res)

	if len(res) != len(exp) {
		t.Error(fmt.Errorf("Number of different mutations is not the same between results and expected %d vs. %d", len(res), len(exp)))
	}

	for i, e := range exp {
		if e != res[i] {
			t.Error(fmt.Errorf("Mutation is not the same between results and expected %s vs. %s", res[i], e))
		}
	}

	if len(resdiffs) != len(expdiffs) {
		t.Error(fmt.Errorf("Number of different mutations is not the same between results and expected %d vs. %d", len(res), len(exp)))
	}

	for i, seqdiff := range expdiffs {
		if len(seqdiff) != len(resdiffs[i]) {
			t.Error(fmt.Errorf("Number of different mutations is not the same between results and expected for seq %d : %d vs. %d", i, len(resdiffs[i]), len(seqdiff)))
		}
		for k, v := range seqdiff {
			if resdiffs[i][k] != v {
				t.Error(fmt.Errorf("Mutation %s does not have the same nb of occurences between results and expected in seq %d : %d vs. %d", k, i, resdiffs[i][k], v))
			}
		}
	}

}

func TestCompress(t *testing.T) {
	in := NewAlign(UNKNOWN)
	in.AddSequence("Seq0000", "GGTTTTTTTT", "")
	in.AddSequence("Seq0001", "TTCCCCACCC", "")
	in.AddSequence("Seq0002", "GGTTTTTTTT", "")
	in.AddSequence("Seq0003", "GGTTTTTTTT", "")
	in.AutoAlphabet()

	exp := NewAlign(UNKNOWN)
	exp.AddSequence("Seq0000", "GTT", "")
	exp.AddSequence("Seq0001", "TAC", "")
	exp.AddSequence("Seq0002", "GTT", "")
	exp.AddSequence("Seq0003", "GTT", "")
	exp.AutoAlphabet()

	expw := []int{2, 1, 7}

	w := in.Compress()

	if len(w) != len(expw) {
		t.Error(fmt.Errorf("Number patterns is not what is expected %v vs. %v", w, expw))
	}

	for i, e := range w {
		if e != expw[i] {
			t.Error(fmt.Errorf("Pattern %d is different %v vs. %v", i, w, expw))
		}
	}

	if !exp.Identical(in) {
		t.Error(fmt.Errorf("Compressed alignment is different from expected \n %s \n vs. \n %s", in.String(), exp.String()))
	}
}

func TestCompress2(t *testing.T) {
	in := NewAlign(UNKNOWN)
	in.AddSequence("Seq0000", "GGGGGGGGGGGGGGGGGGGG", "")
	in.AddSequence("Seq0001", "TTTTTTTTTTTTTTTTTTTT", "")
	in.AddSequence("Seq0002", "GGGGGGGGGGGGGGGGGGGG", "")
	in.AddSequence("Seq0003", "AAAAAAAAAAAAAAAAAAAA", "")
	in.AutoAlphabet()

	exp := NewAlign(UNKNOWN)
	exp.AddSequence("Seq0000", "G", "")
	exp.AddSequence("Seq0001", "T", "")
	exp.AddSequence("Seq0002", "G", "")
	exp.AddSequence("Seq0003", "A", "")
	exp.AutoAlphabet()

	expw := []int{20}

	w := in.Compress()

	if len(w) != len(expw) {
		t.Error(fmt.Errorf("Number patterns is not what is expected %v vs. %v", w, expw))
	}

	for i, e := range w {
		if e != expw[i] {
			t.Error(fmt.Errorf("Pattern %d is different %v vs. %v", i, w, expw))
		}
	}

	if !exp.Identical(in) {
		t.Error(fmt.Errorf("Compressed alignment is different from expected \n %s \n vs. \n %s", in.String(), exp.String()))
	}
}
