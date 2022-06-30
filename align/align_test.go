package align

import (
	"fmt"
	"math"
	"reflect"
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

	a.IterateChar(func(name string, sequence []uint8) bool {
		if !strings.HasPrefix(name, "IDENT") {
			t.Error("Sequence name does not start with expected id: IDENT")
			return true
		}
		return false
	})

	a.AppendSeqIdentifier("IDENT", true)
	a.IterateChar(func(name string, sequence []uint8) bool {
		if !strings.HasSuffix(name, "IDENT") {
			t.Error("Sequence name does not end with expected id: IDENT")
			return true
		}
		return false
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
	a.IterateChar(func(name string, sequence []uint8) bool {
		expected, found := a2.GetSequenceNameById(i)
		if !found {
			t.Error("Unknown sequence name after clean")
			return true
		}
		if name != expected {
			t.Error("Unexpected sequence name after clean")
			return true
		}
		i++
		return false
	})
}

func TestRemoveOneGapSite(t *testing.T) {
	var start, end int
	var kept []int
	var l int
	a, err := RandomAlignment(AMINOACIDS, 300, 300)
	if err != nil {
		t.Error(err)

	}

	l = a.Length()
	/* We add 1 gap per site */
	pos := 0
	a.IterateChar(func(name string, sequence []uint8) bool {
		sequence[pos] = GAP
		pos++
		return false
	})

	start, end, kept = a.RemoveGapSites(0.0, false)

	if start != l {
		t.Errorf("We should have removed all positions from start: %d %d", start, l)
	}

	if end != l {
		t.Errorf("We should have removed all positions from end: %d %d", end, l)
	}

	if a.Length() != 0 {
		t.Error("We should have removed all positions")
	}

	if len(kept) != 0 {
		t.Error("We should have removed all positions in kept slice")
	}

	a.IterateChar(func(name string, sequence []uint8) bool {
		if len(sequence) != 0 {
			t.Error(fmt.Sprintf("Sequence length after removing gaps should be 0 and is : %d", len(sequence)))
			return true
		}
		return false
	})
}

func TestRemoveAllGapSites(t *testing.T) {
	var start, end int
	var kept []int

	a, err := RandomAlignment(AMINOACIDS, 300, 300)
	if err != nil {
		t.Error(err)

	}

	backupseq := make([]uint8, 0, 300)
	seq0, found := a.GetSequenceCharById(0)
	if !found {
		t.Error("Problem finding first sequence")
	}

	/* We add all gaps on 1 site */
	/* And one gap at all sites */
	pos1 := 20
	pos2 := 0
	a.IterateChar(func(name string, sequence []uint8) bool {
		sequence[pos1] = GAP
		sequence[pos2] = GAP
		pos2++
		return false
	})
	backupseq = append(backupseq, seq0...)
	/* Remove position 20 */
	backupseq = append(backupseq[:20], backupseq[21:]...)

	start, end, kept = a.RemoveGapSites(1.0, false)

	if start != 0 {
		t.Errorf("We should have removed 0 positions from start: %d", start)
	}

	if end != 0 {
		t.Errorf("We should have removed 0 positions from end: %d", end)
	}

	if a.Length() != 299 {
		t.Error("We should have removed only one position")
	}

	if len(kept) != 299 {
		t.Errorf("We should have removed only one position in kept slice (%d remaining)", len(kept))
	}

	a.IterateChar(func(name string, sequence []uint8) bool {
		if len(sequence) != 299 {
			t.Error(fmt.Sprintf("Sequence length after removing gaps should be 299 and is : %d", len(sequence)))
			return true
		}
		return false
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

func TestRemoveOneGapSiteEnds(t *testing.T) {
	var l, start, end int
	var kept []int

	a, err := RandomAlignment(AMINOACIDS, 300, 300)
	if err != nil {
		t.Error(err)

	}

	l = a.Length()
	/* We add 1 gap per site */
	pos := 0
	a.IterateChar(func(name string, sequence []uint8) bool {
		sequence[pos] = GAP
		pos++
		return false
	})

	start, end, kept = a.RemoveGapSites(0.0, true)

	if start != l {
		t.Errorf("We should have removed all positions from start: %d %d", start, l)
	}

	if end != l {
		t.Errorf("We should have removed all positions from end: %d %d", end, l)
	}

	if a.Length() != 0 {
		t.Error("We should have removed all positions")
	}

	if len(kept) != 0 {
		t.Error("We should have removed all positions in kept slice")
	}

	a.IterateChar(func(name string, sequence []uint8) bool {
		if len(sequence) != 0 {
			t.Error(fmt.Sprintf("Sequence length after removing gaps should be 0 and is : %d", len(sequence)))
		}
		return false
	})
}

func TestRemoveAllGapSitesEnds(t *testing.T) {
	var start, end int
	var kept []int

	a, err := RandomAlignment(AMINOACIDS, 300, 300)
	if err != nil {
		t.Error(err)

	}

	/* We add all gaps on 1 site */
	/* And one gap at all sites */
	pos1 := 20
	pos2 := 0
	a.IterateChar(func(name string, sequence []uint8) bool {
		sequence[pos1] = GAP
		sequence[pos2] = GAP
		pos2++
		return false
	})

	start, end, kept = a.RemoveGapSites(1.0, true)

	if start != 0 {
		t.Errorf("We should have removed 0 positions from start: %d", start)
	}

	if end != 0 {
		t.Errorf("We should have removed 0 positions from end: %d", end)
	}

	if a.Length() != 300 {
		t.Error(fmt.Sprintf("We should not have removed any positions: %d", a.Length()))
	}

	if len(kept) != 300 {
		t.Error("We should not have removed any position in kept slice")
	}

	a.IterateChar(func(name string, sequence []uint8) bool {
		if len(sequence) != 300 {
			t.Error(fmt.Sprintf("Sequence length after removing gaps should be 299 and is : %d", len(sequence)))
			return true
		}
		return false
	})
}

func TestRemoveGapSitesEnds(t *testing.T) {
	var start, end int
	var kept []int

	in := NewAlign(UNKNOWN)
	in.AddSequence("Seq0000", "--GGTCCACTCTTTTGTCTT-TACCTA-G--", "")
	in.AddSequence("Seq0001", "G---CACCGGC-CGTAATGACG-ACCC--T-", "")
	in.AddSequence("Seq0002", "-T-G-TTTCCTGC-AACAT-ACC-AAC-C--", "")
	in.AddSequence("Seq0003", "A-ACCACAACAGTCA-GTACTCTT-TG--T-", "")
	in.AddSequence("Seq0004", "-----GAAGG-CCAAGGT-TCGCCGCCC---", "")
	in.AutoAlphabet()

	exp := NewAlign(UNKNOWN)
	exp.AddSequence("Seq0000", "GTCCACTCTTTTGTCTT-TACCTA", "")
	exp.AddSequence("Seq0001", "-CACCGGC-CGTAATGACG-ACCC", "")
	exp.AddSequence("Seq0002", "G-TTTCCTGC-AACAT-ACC-AAC", "")
	exp.AddSequence("Seq0003", "CCACAACAGTCA-GTACTCTT-TG", "")
	exp.AddSequence("Seq0004", "--GAAGG-CCAAGGT-TCGCCGCC", "")
	exp.AutoAlphabet()

	start, end, kept = in.RemoveGapSites(0.5, true)

	if start != 3 {
		t.Errorf("We should have removed 3 positions from start: %d", start)
	}

	if end != 4 {
		t.Errorf("We should have removed 4 positions from end: %d", end)
	}

	if len(kept) != 24 {
		t.Errorf("We should have kept 24 positions in kept slice, not %d", len(kept))
	}

	if !exp.Identical(in) {
		t.Error(fmt.Errorf("Alignment after removing gap ends is different from expected \n %s \n vs. \n %s", in.String(), exp.String()))
	}
}

func TestRemoveOneGapSequence(t *testing.T) {
	a, err := RandomAlignment(AMINOACIDS, 300, 300)
	if err != nil {
		t.Error(err)

	}

	/* We add 1 gap per site on all sequences*/
	pos := 0
	a.IterateChar(func(name string, sequence []uint8) bool {
		sequence[pos] = GAP
		pos++
		return false
	})

	a.RemoveGapSeqs(0.0, false)

	if a.NbSequences() != 0 {
		t.Error(fmt.Errorf("We should have removed all sequences %d", a.NbSequences()))
	}
}

func TestRemoveOneGapSequence2(t *testing.T) {
	a, err := RandomAlignment(AMINOACIDS, 300, 300)
	if err != nil {
		t.Error(err)

	}

	/* We add 1 gap per site on half of the sequences*/
	pos := 0
	a.IterateChar(func(name string, sequence []uint8) bool {
		if pos%2 == 0 {
			sequence[pos] = GAP
		}
		pos++
		return false
	})

	a.RemoveGapSeqs(0.0, false)

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

	a.RemoveGapSeqs(1.0, false)

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

	a.RemoveGapSeqs(0.5, false)

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
	a.IterateChar(func(name string, sequence []uint8) bool {
		sequence[pos] = GAP
		pos++
		return false
	})

	a2, err2 := a.Clone()
	if err2 != nil {
		t.Error(err2)
	}

	a.RemoveGapSites(0.0, false)

	a2.IterateChar(func(name string, sequence []uint8) bool {
		if len(sequence) != 300 {
			t.Error(fmt.Sprintf("Clone lenght should be 300 and is : %d", len(sequence)))
			return true
		}
		return false
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
	a2.IterateChar(func(name string, sequence []uint8) bool {
		s2, ok := a.GetSequenceCharById(i)
		n2, ok2 := a.GetSequenceNameById(i)

		if !ok || !ok2 {
			t.Error(fmt.Sprintf("Sequence not found in clone alignment: %s", name))
			return true
		}

		if len(sequence) != len(s2) {
			t.Error(fmt.Sprintf("Clone length is different from original length : %d != %d", len(sequence), len(s2)))
			return true
		}
		if name != n2 {
			t.Error(fmt.Sprintf("Clone and original sequences at position %d have different names : %s != %s", i, name, n2))
			return true
		}
		for j, c := range sequence {
			if c != s2[j] {
				t.Error(fmt.Sprintf("Clone sequence is different from original at position %d : %c != %c", j, c, s2[j]))
				return true
			}
		}
		i++
		return false
	})
}

func TestAvgAlleles(t *testing.T) {
	a, err := RandomAlignment(AMINOACIDS, 300, 300)
	if err != nil {
		t.Error(err)

	}

	a.IterateChar(func(name string, sequence []uint8) bool {
		for j := range sequence {
			sequence[j] = 'A'
		}
		return false
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
	a.IterateChar(func(name string, sequence []uint8) bool {
		for j := range sequence {
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
		return false
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
	a.IterateChar(func(name string, sequence []uint8) bool {
		expected := fmt.Sprintf("New%04d", i)
		if name != expected {
			t.Error(fmt.Sprintf("Sequence name should be %s and is %s", expected, name))
			return true
		}
		i++
		return false
	})

	a.Rename(namemap2)
	i = 0
	a.IterateChar(func(name string, sequence []uint8) bool {
		expected := fmt.Sprintf("Seq%04d", i)
		if name != expected {
			t.Error(fmt.Sprintf("Sequence name should be %s and is %s", expected, name))
			return true
		}
		i++
		return false
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
	a.IterateChar(func(name string, sequence []uint8) bool {
		for _, c := range sequence {
			if c == 'A' {
				acount++
			}
		}
		return false
	})

	gapcount := 0
	a.Replace("A", "-", false)
	a.IterateChar(func(name string, sequence []uint8) bool {
		for _, c := range sequence {
			if c == 'A' {
				t.Errorf("There should not remains A after replace")
				return true
			}
			if c == '-' {
				gapcount++
			}
		}
		return false
	})
	if gapcount != acount {
		t.Errorf("Each A should have been replaced by a -")
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
	a.IterateChar(func(name string, sequence []uint8) bool {
		prev := uint8('0')
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
		return false
	})

	gapcount := 0
	a.Replace("A.", "--", true)
	a.IterateChar(func(name string, sequence []uint8) bool {
		prev := uint8('0')
		for _, c := range sequence {
			if (c == 'A' && prev == 'A') ||
				(c == 'C' && prev == 'A') ||
				(c == 'G' && prev == 'A') ||
				(c == 'T' && prev == 'A') {
				t.Error(fmt.Sprintf("There should not remains %c%c after replace", prev, c))
				return true
			}
			if c == '-' {
				gapcount++
			}
			prev = c
		}
		return false
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
	a, _ := RandomAlignment(AMINOACIDS, length, nbseqs)

	alldifferent := []uint8{'A', 'R', 'N', 'D', 'C'}
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
	a, _ := RandomAlignment(AMINOACIDS, length, nbseqs)

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
	subalign.IterateChar(func(name string, sequence []uint8) bool {
		if len(sequence) != 90 {
			t.Error(fmt.Sprintf("Length of subsequence must be %d and is %d", 90, len(sequence)))
			return true
		}
		for i, a := range sequence {
			if a != 'A' {
				t.Error(fmt.Sprintf("Character at position %d must be %c and is %c", i, 'A', a))
				return true
			}
		}
		return false
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

	a.IterateChar(func(name string, sequence []uint8) bool {
		s, _ := acopy.GetSequence(name)
		s2, _ := a2.GetSequence(name)
		if string(sequence) != s+s2 {
			t.Errorf("Concatenated sequence is not correct")
			return true
		}
		return false
	})
}

func TestAppend(t *testing.T) {
	var err error
	var a, a2, a3 Alignment

	a = NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "ACGT", "")
	a.AddSequence("B", "ACGG", "")
	a.AddSequence("C", "ACGT", "")
	a.AddSequence("C", "ACGT", "")
	a.AddSequence("D", "ACGT", "")

	a2 = NewAlign(NUCLEOTIDS)
	a2.AddSequence("E", "ACGT", "")
	a2.AddSequence("F", "ACGG", "")
	a2.AddSequence("G", "ACGT", "")
	a2.AddSequence("H", "ACGT", "")
	a2.AddSequence("I", "ACGT", "")

	a3 = NewAlign(NUCLEOTIDS)
	a3.AddSequence("A", "ACGT", "")
	a3.AddSequence("B", "ACGG", "")
	a3.AddSequence("C", "ACGT", "")
	a3.AddSequence("C", "ACGT", "")
	a3.AddSequence("D", "ACGT", "")
	a3.AddSequence("E", "ACGT", "")
	a3.AddSequence("F", "ACGG", "")
	a3.AddSequence("G", "ACGT", "")
	a3.AddSequence("H", "ACGT", "")
	a3.AddSequence("I", "ACGT", "")

	if err = a.Append(a2); err != nil {
		t.Error(err)
	}

	if !a.Identical(a3) {
		t.Errorf("After append, alignment is not expected")
	}
}
func TestDedup(t *testing.T) {
	var err error
	var a Alignment = NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "ACGT", "")
	a.AddSequence("B", "ACGG", "")
	a.AddSequence("C", "ACGT", "")
	a.AddSequence("C", "ACGT", "")
	a.AddSequence("D", "ACGT", "")

	var expectedIdentical [][]string = [][]string{{"A", "C", "C_0001", "D"}, {"B"}}
	var identical [][]string

	if a.NbSequences() != 5 {
		t.Error("There should be 5 sequences before deduplicaion")
	}

	if identical, err = a.Deduplicate(); err != nil {
		t.Error(err)
	}

	if a.NbSequences() != 2 {
		t.Error("There should be 2 sequences in the deduplicated alignment")
	}

	if len(identical) != len(expectedIdentical) {
		t.Errorf("After deduplicate, slice of identical sequences is different from expected: %d vs. %d", len(expectedIdentical), len(identical))
		return
	}
	for i, id := range expectedIdentical {
		if len(id) != len(identical[i]) {
			t.Errorf("After deduplicate, slice of identical sequences is different from expected: %d vs. %d", len(id), len(identical[i]))
			return
		}
		for j, name := range id {
			if identical[i][j] != name {
				t.Errorf("After deduplicate, slice of identical sequences is different from expected: %s vs. %s", name, identical[i][j])
			}
		}
	}

	if _, err = a.Deduplicate(); err != nil {
		t.Error(err)
	}
	if a.NbSequences() != 2 {
		t.Error("There should be 2 sequences in the de-deduplicated alignment")
	}
}

func TestDedup2(t *testing.T) {
	var err error
	var sb SeqBag = NewSeqBag(NUCLEOTIDS)
	sb.AddSequence("A", "ACGT", "")
	sb.AddSequence("B", "ACGG", "")
	sb.AddSequence("C", "ACGT", "")
	sb.AddSequence("C", "ACGT", "")
	sb.AddSequence("D", "ACGT", "")

	var expectedIdentical [][]string = [][]string{{"A", "C", "C_0001", "D"}, {"B"}}
	var identical [][]string

	if sb.NbSequences() != 5 {
		t.Error("There should be 5 sequences before deduplicaion")
	}

	if identical, err = sb.Deduplicate(); err != nil {
		t.Error(err)
	}

	if sb.NbSequences() != 2 {
		t.Error("There should be 2 sequences in the deduplicated alignment")
	}

	if len(identical) != len(expectedIdentical) {
		t.Errorf("After deduplicate, slice of identical sequences is different from expected: %d vs. %d", len(expectedIdentical), len(identical))
		return
	}
	for i, id := range expectedIdentical {
		if len(id) != len(identical[i]) {
			t.Errorf("After deduplicate, slice of identical sequences is different from expected: %d vs. %d", len(id), len(identical[i]))
			return
		}
		for j, name := range id {
			if identical[i][j] != name {
				t.Errorf("After deduplicate, slice of identical sequences is different from expected: %s vs. %s", name, identical[i][j])
			}
		}
	}

	if _, err = sb.Deduplicate(); err != nil {
		t.Error(err)
	}
	if sb.NbSequences() != 2 {
		t.Error("There should be 2 sequences in the de-deduplicated alignment")
	}
}

func TestDedup3(t *testing.T) {
	var err error
	var sb SeqBag = NewSeqBag(NUCLEOTIDS)
	sb.AddSequence("A", "ACGT", "")
	sb.AddSequence("B", "ACGG", "")
	sb.AddSequence("C", "ACGT", "")
	sb.AddSequence("C", "ACG", "")
	sb.AddSequence("D", "A", "")

	var expectedIdentical [][]string = [][]string{{"A", "C"}, {"B"}, {"C_0001"}, {"D"}}
	var identical [][]string

	if sb.NbSequences() != 5 {
		t.Error("There should be 5 sequences before deduplicaion")
	}

	if identical, err = sb.Deduplicate(); err != nil {
		t.Error(err)
	}

	if sb.NbSequences() != 4 {
		t.Error("There should be 4 sequences in the deduplicated alignment")
	}

	if len(identical) != len(expectedIdentical) {
		t.Errorf("After deduplicate, slice of identical sequences is different from expected: %d vs. %d", len(expectedIdentical), len(identical))
		return
	}
	for i, id := range expectedIdentical {
		if len(id) != len(identical[i]) {
			t.Errorf("After deduplicate, slice of identical sequences is different from expected: %d vs. %d", len(id), len(identical[i]))
			return
		}
		for j, name := range id {
			if identical[i][j] != name {
				t.Errorf("After deduplicate, slice of identical sequences is different from expected: %s vs. %s", name, identical[i][j])
			}
		}
	}

	if _, err = sb.Deduplicate(); err != nil {
		t.Error(err)
	}
	if sb.NbSequences() != 4 {
		t.Error("There should be 4 sequences in the de-deduplicated alignment")
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

func TestTranslateMitoI(t *testing.T) {

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
	exp3phases.AddSequence("Seq0000_1", "WMTFSVALPHRNSLPHRVAKSLPQRVAEGGITITLS*LS*CNDS", "")
	exp3phases.AddSequence("Seq0000_2", "GWLFLLPSPTATVFLIE*RSDYHNGWRSVASLLHYHSCRSVMM", "")
	exp3phases.AddSequence("Seq0001_0", "MDDFFCCPPPPQQSSSSSSEETTTTGGGGWHHYYIIMVVV*W", "")
	exp3phases.AddSequence("Seq0001_1", "WMTFSVALPHRNSLPHRVAKSLPQRVAEGGITITLS*LSYN", "")
	exp3phases.AddSequence("Seq0001_2", "GWLFLLPSPTATVFLIE*RSDYHNGWRSVASLLHYHSCRMM", "")
	exp3phases.AutoAlphabet()

	if err := in.Translate(0, 2); err != nil {
		t.Error(err)
	} else {
		if !expaa.Identical(in) {
			t.Error(fmt.Errorf("Expected sequences are different from phased sequences"))
		}
	}

	if err := in2.Translate(-1, 2); err != nil {
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

	exp2 := NewSeqBag(UNKNOWN)
	exp2.AddSequence("Seq0000", "GATTAATT--CCGTAGGCCA", "")
	exp2.AddSequence("Seq0001", "GAATCTGA--ATCGAACACT", "")
	exp2.AddSequence("Seq0002", "TTAAGTTT--ACTTCTAATG", "")
	exp2.AutoAlphabet()

	exp3 := NewSeqBag(UNKNOWN)
	exp3.AddSequence("Seq0000", "GATTAATT,,CCGTAGGCCA", "")
	exp3.AddSequence("Seq0001", "GAATCTGA,,ATCGAACACT", "")
	exp3.AddSequence("Seq0002", "TTAAGTTT,,ACTTCTAATG", "")
	exp3.AutoAlphabet()

	exp4 := NewAlign(UNKNOWN)
	exp4.AddSequence("Seq0000", "GATTAATTTGACGTAGGCCA", "")
	exp4.AddSequence("Seq0001", "GAATCTGAAGACCGAACACT", "")
	exp4.AddSequence("Seq0002", "TTAAGTTTTGACTTCTAATG", "")
	exp4.AutoAlphabet()

	res, _ := in.Clone()
	if err := res.Mask(8, 2, "", false); err != nil {
		t.Error(err)
	} else {
		if !exp.Identical(res) {
			t.Error(fmt.Errorf("Expected sequences are different from masked sequences"))
		}
	}

	res, _ = in.Clone()
	if err := res.Mask(8, 2, "AMBIG", false); err != nil {
		t.Error(err)
	} else {
		if !exp.Identical(res) {
			t.Error(fmt.Errorf("Expected sequences are different from masked sequences"))
		}
	}

	res, _ = in.Clone()
	if err := res.Mask(8, 2, "GAP", false); err != nil {
		t.Error(err)
	} else {
		if !exp2.Identical(res) {
			t.Error(fmt.Errorf("Expected sequences are different from masked sequences"))
		}
	}

	res, _ = in.Clone()
	if err := res.Mask(8, 2, ",", false); err != nil {
		t.Error(err)
	} else {
		if !exp3.Identical(res) {
			t.Error(fmt.Errorf("Expected sequences are different from masked sequences"))
		}
	}

	res, _ = in.Clone()
	if err := res.Mask(9, 3, "MAJ", false); err != nil {
		t.Error(err)
	} else {
		if !exp4.Identical(res) {
			fmt.Println(exp4)
			fmt.Println(res)
			t.Error(fmt.Errorf("Expected sequences are different from masked sequences"))
		}
	}
}

func TestMaskNoGAP(t *testing.T) {
	in := NewAlign(UNKNOWN)
	in.AddSequence("Seq0000", "GATTAATTTGCCGTAGGCCA", "")
	in.AddSequence("Seq0001", "GAATCTGAA-ATCGAACACT", "")
	in.AddSequence("Seq0002", "TTAAGTTTT-ACTTCTAATG", "")
	in.AutoAlphabet()

	exp := NewSeqBag(UNKNOWN)
	exp.AddSequence("Seq0000", "GATTAATTTNCCGTAGGCCA", "")
	exp.AddSequence("Seq0001", "GAATCTGAANATCGAACACT", "")
	exp.AddSequence("Seq0002", "TTAAGTTTTNACTTCTAATG", "")
	exp.AutoAlphabet()

	exp2 := NewSeqBag(UNKNOWN)
	exp2.AddSequence("Seq0000", "GATTAATTT-CCGTAGGCCA", "")
	exp2.AddSequence("Seq0001", "GAATCTGAA-ATCGAACACT", "")
	exp2.AddSequence("Seq0002", "TTAAGTTTT-ACTTCTAATG", "")
	exp2.AutoAlphabet()

	exp3 := NewSeqBag(UNKNOWN)
	exp3.AddSequence("Seq0000", "GATTAATTTNCCGTAGGCCA", "")
	exp3.AddSequence("Seq0001", "GAATCTGAA-ATCGAACACT", "")
	exp3.AddSequence("Seq0002", "TTAAGTTTT-ACTTCTAATG", "")
	exp3.AutoAlphabet()

	res, _ := in.Clone()
	if err := res.Mask(9, 1, "", false); err != nil {
		t.Error(err)
	} else {
		if !exp.Identical(res) {
			t.Error(fmt.Errorf("Expected sequences are different from masked sequences"))
		}
	}

	res, _ = in.Clone()
	if err := res.Mask(9, 1, "MAJ", true); err != nil {
		t.Error(err)
	} else {
		if !exp2.Identical(res) {
			t.Error(fmt.Errorf("Expected sequences are different from masked sequences"))
		}
	}

	res, _ = in.Clone()
	if err := res.Mask(9, 1, "", true); err != nil {
		t.Error(err)
	} else {
		if !exp3.Identical(res) {
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

	if err := in.Mask(12, 2000, "AMBIG", false); err != nil {
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

func TestDiffReverse(t *testing.T) {
	in := NewAlign(UNKNOWN)
	in.AddSequence("Seq0000", "GATTAATTTGCCGTAGGCCAGAATCTGAAGATCGAACACTTTAAGTTTTCACTTCTAATGGAGAGGACTAGTTCATACTTTTTAAACACTTTTACATCGA", "")
	in.AddSequence("Seq0001", "TG.CGGACCTAA...TTGAGT.CAAC.GT.TATTCCAG.GG.GGAGAGGTCTA.T.TTCC.GTT.A.GG.C.CT.G.GC.G.A..GGGTA.GGC...GTG", "")
	in.AddSequence("Seq0002", "CTAAGCGCG.G..G.TTG.T.TTGGA.C.AGGTT..ATAC.CGGCAA.G.C.CATG.TCCCCC.A.GAC.A.AAGAG.GAAG.T.GA..AAA.GA.C.CC", "")
	in.AddSequence("Seq0003", "..G.GGAGGCTTTAT...ACA.GGTATT...GACTGAGGGGC.CCCCGG..TGGTA.GCA.GAGCC.TCGCGAAGGCT.CAGGT.T.TTCC.GTGT.ACC", "")
	in.AutoAlphabet()

	out, _ := in.Clone()
	out.ReplaceMatchChars()

	i := 0
	var ref []uint8
	out.IterateChar(func(name string, seq []uint8) bool {
		orig, _ := in.GetSequenceCharById(i)
		if i == 0 {
			ref = seq
			for site := 0; site < out.Length(); site++ {
				if seq[site] != orig[site] {
					t.Errorf("Original reference character has been changed (%c vs. %c)", orig[site], seq[site])
					return true
				}
			}
		} else {
			for site := 0; site < out.Length(); site++ {
				if orig[site] == POINT {
					if seq[site] != ref[site] {
						t.Errorf(". has not been replaced by the right character (%c vs. %c)", ref[site], seq[site])
						return true
					}
				} else {
					if seq[site] != orig[site] {
						t.Errorf("Original character has been changed (%c vs. %c)", orig[site], seq[site])
						return true
					}
				}
			}
		}
		i++
		return false
	})
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

func Test_align_RefCoordinates(t *testing.T) {
	var err error
	var alistart, alilen int

	in := NewAlign(UNKNOWN)

	in.AddSequence("Seq0000", "--ACG--AT---GC", "")
	in.AddSequence("Seq0001", "GGACGTTATCGGGC", "")
	in.AutoAlphabet()

	expstart, explen := 2, 6
	alistart, alilen, err = in.RefCoordinates("Seq0000", 0, 4)
	if alistart != expstart || alilen != explen || err != nil {
		t.Error(fmt.Errorf("alistart: %d!=%d | alilen: %d!=%d | err: %v", alistart, expstart, alilen, explen, err))
	}

	expstart, explen = 3, 6
	alistart, alilen, err = in.RefCoordinates("Seq0000", 1, 4)
	if alistart != expstart || alilen != explen || err != nil {
		t.Error(fmt.Errorf("alistart: %d!=%d | alilen: %d!=%d | err: %v", alistart, expstart, alilen, explen, err))
	}

	expstart, explen = 3, 11
	alistart, alilen, err = in.RefCoordinates("Seq0000", 1, 6)
	if alistart != expstart || alilen != explen || err != nil {
		t.Error(fmt.Errorf("alistart: %d!=%d | alilen: %d!=%d | err: %v", alistart, expstart, alilen, explen, err))
	}

	expstart, explen = 1, 6
	alistart, alilen, err = in.RefCoordinates("Seq0001", 1, 6)
	if alistart != expstart || alilen != explen || err != nil {
		t.Error(fmt.Errorf("alistart: %d!=%d | alilen: %d!=%d | err: %v", alistart, expstart, alilen, explen, err))
	}

	expstart, explen = 12, 2
	alistart, alilen, err = in.RefCoordinates("Seq0000", 5, 2)
	if alistart != expstart || alilen != explen || err != nil {
		t.Error(fmt.Errorf("alistart: %d!=%d | alilen: %d!=%d | err: %v", alistart, expstart, alilen, explen, err))
	}

}

func TestConsensus(t *testing.T) {
	var a Alignment
	var c Alignment
	var exp Alignment

	a = NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "ACGACGACGACC", "")
	a.AddSequence("B", "ATCTT-TTTTTC", "")
	a.AddSequence("C", "ATCTT-TTTTTT", "")

	exp = NewAlign(NUCLEOTIDS)
	exp.AddSequence("consensus", "ATCTT-TTTTTC", "")

	c = a.Consensus(false, false)

	if !exp.Identical(c) {
		t.Error(fmt.Errorf("Consensus is not identical to expected alignment"))
	}
}

func TestConsensusGaps(t *testing.T) {
	var a Alignment
	var c Alignment
	var exp Alignment

	a = NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "ACGACGACGACC", "")
	a.AddSequence("B", "ATCTT-TTTTTC", "")
	a.AddSequence("C", "ATCTT-TTTTTT", "")

	exp = NewAlign(NUCLEOTIDS)
	exp.AddSequence("consensus", "ATCTTGTTTTTC", "")

	c = a.Consensus(true, false)

	if !exp.Identical(c) {
		fmt.Println(exp)
		fmt.Println(c)
		t.Error(fmt.Errorf("Consensus is not identical to expected alignment"))
	}
}

func TestConsensusNs(t *testing.T) {
	var a Alignment
	var c Alignment
	var exp Alignment

	a = NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "ACGACGACGACCN", "")
	a.AddSequence("B", "ATCTTNTTTTT-N", "")
	a.AddSequence("C", "ATCTTNTTTTT-T", "")

	exp = NewAlign(NUCLEOTIDS)
	exp.AddSequence("consensus", "ATCTTGTTTTT-T", "")

	c = a.Consensus(false, true)

	if !exp.Identical(c) {
		fmt.Println(exp)
		fmt.Println(c)
		t.Error(fmt.Errorf("Consensus is not identical to expected alignment"))
	}
}

func Test_align_NumGapsUniquePerSequence(t *testing.T) {
	var a, ref Alignment
	var ng, nn, nb []int
	var expng, expnn, expnb []int
	var err error
	var profile *CountProfile

	a = NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "ACGACGA-GACC", "")
	a.AddSequence("B", "AT-TT-T-TTTC", "")
	a.AddSequence("C", "ATCTT-TTT--T", "")

	ref = NewAlign(NUCLEOTIDS)
	ref.AddSequence("A", "AAAAAAAAAAAA", "")
	ref.AddSequence("B", "AAAAAAAAAAAA", "")
	ref.AddSequence("C", "AAAAAAAAAAAA", "")

	expng = []int{0, 1, 2}
	expnn = []int{1, 3, 3}
	expnb = []int{0, 1, 2}

	profile = NewCountProfileFromAlignment(ref)
	if ng, nn, nb, err = a.NumGapsUniquePerSequence(profile); err != nil {
		t.Error(err)
	}

	if !reflect.DeepEqual(expng, ng) {
		t.Error(fmt.Errorf("Numgaps is not what is expected, have %v, want %v", ng, expng))
	}
	if !reflect.DeepEqual(expnn, nn) {
		t.Error(fmt.Errorf("Numgaps is not what is expected, have %v, want %v", nn, expnn))
	}
	if !reflect.DeepEqual(expnb, nb) {
		t.Error(fmt.Errorf("Numgaps is not what is expected, have %v, want %v", nb, expnb))
	}
}

func Test_align_NumMutationsUniquePerSequence(t *testing.T) {
	var a, ref Alignment
	var ng, nn, nb []int
	var expng, expnn, expnb []int
	var err error
	var profile *CountProfile

	a = NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "ACGACGA-GACC", "")
	a.AddSequence("B", "AT-TT-T-TTTC", "")
	a.AddSequence("C", "ATCTT-TTT--T", "")

	ref = NewAlign(NUCLEOTIDS)
	ref.AddSequence("A", "AAAAAAAAAAAA", "")
	ref.AddSequence("B", "AAAAAAAAAAAA", "")
	ref.AddSequence("C", "AAAAAAAAAAAA", "")

	expng = []int{9, 2, 3}
	expnn = []int{7, 8, 8}
	expnb = []int{6, 2, 3}

	profile = NewCountProfileFromAlignment(ref)
	if ng, nn, nb, err = a.NumMutationsUniquePerSequence(profile); err != nil {
		t.Error(err)
	}

	if !reflect.DeepEqual(expng, ng) {
		t.Error(fmt.Errorf("Numgaps is not what is expected, have %v, want %v", ng, expng))
	}
	if !reflect.DeepEqual(expnn, nn) {
		t.Error(fmt.Errorf("Numgaps is not what is expected, have %v, want %v", nn, expnn))
	}

	if !reflect.DeepEqual(expnb, nb) {
		t.Error(fmt.Errorf("Numgaps is not what is expected, have %v, want %v", nb, expnb))
	}
}

func Test_align_NumMutationsComparedToReferenceSequence(t *testing.T) {
	var a Alignment
	var ng int
	var exp []int
	var err error

	a = NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "ACGACGA-GACC", "")
	a.AddSequence("B", "AT-TT-T-TTTC", "")
	a.AddSequence("C", "ATCTT-TTT--T", "")
	a.AddSequence("D", "CCCCCCCCCCCC", "")

	ref := NewSequence("ref", []uint8("CCCCCCCCCCCC"), "")

	exp = []int{7, 8, 8, 0}

	for i, s := range a.Sequences() {
		if ng, err = s.NumMutationsComparedToReferenceSequence(a.Alphabet(), ref); err != nil {
			t.Error(err)
		}
		if !reflect.DeepEqual(exp[i], ng) {
			t.Error(fmt.Errorf("Nummutations is not what is expected, have %v, want %v", ng, exp))
		}
	}
}

func Test_align_RemoveMajorityCharacterSites(t *testing.T) {
	var a, a2, a3, a4, exp, exp2, exp3, exp4 Alignment

	a = NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "N-ANGA-GACNC", "")
	a.AddSequence("B", "N-TN-T-TTTAC", "")
	a.AddSequence("C", "NCTN-TTT--AT", "")
	a.AddSequence("D", "C-ANCCCCCCCC", "")

	a2, _ = a.Clone()
	a3, _ = a.Clone()
	a4, _ = a.Clone()

	exp = NewAlign(NUCLEOTIDS)
	exp.AddSequence("A", "AGA-GACN", "")
	exp.AddSequence("B", "T-T-TTTA", "")
	exp.AddSequence("C", "T-TTT--A", "")
	exp.AddSequence("D", "ACCCCCCC", "")

	exp2 = NewAlign(NUCLEOTIDS)
	exp2.AddSequence("A", "ANGA-GACN", "")
	exp2.AddSequence("B", "TN-T-TTTA", "")
	exp2.AddSequence("C", "TN-TTT--A", "")
	exp2.AddSequence("D", "ANCCCCCCC", "")

	exp3 = NewAlign(NUCLEOTIDS)
	exp3.AddSequence("A", "AGA-GAN", "")
	exp3.AddSequence("B", "T-T-TTA", "")
	exp3.AddSequence("C", "T-TTT-A", "")
	exp3.AddSequence("D", "ACCCCCC", "")

	exp4 = NewAlign(NUCLEOTIDS)
	exp4.AddSequence("A", "AGA-GAC", "")
	exp4.AddSequence("B", "T-T-TTT", "")
	exp4.AddSequence("C", "T-TTT--", "")
	exp4.AddSequence("D", "ACCCCCC", "")

	a.RemoveMajorityCharacterSites(0.6, false, false, false)
	a2.RemoveMajorityCharacterSites(0.6, true, false, false)
	a3.RemoveMajorityCharacterSites(0.6, false, true, false)
	a4.RemoveMajorityCharacterSites(0.6, false, false, true)

	if !a.Identical(exp) {
		t.Errorf(a.String())
		t.Error(fmt.Errorf("Remove majority failed"))
	}
	if !a2.Identical(exp2) {
		t.Errorf(a2.String())
		t.Error(fmt.Errorf("Remove majority failed with ends"))
	}
	if !a3.Identical(exp3) {
		t.Errorf(a3.String())
		t.Error(fmt.Errorf("Remove majority failed ignore gaps"))
	}
	if !a4.Identical(exp4) {
		t.Errorf(a4.String())
		t.Error(fmt.Errorf("Remove majority failed ignore Ns"))
	}
}

func Test_align_MaskUnique(t *testing.T) {
	var a, res, exp, exp2, exp3 Alignment

	a = NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "ACANGA-TACC", "")
	a.AddSequence("B", "ACTN-T-TTTC", "")
	a.AddSequence("C", "ACTN-TTT--T", "")
	a.AddSequence("D", "C-ANCCCCCCC", "")

	exp = NewAlign(NUCLEOTIDS)
	exp.AddSequence("A", "ACANNN-TNCC", "")
	exp.AddSequence("B", "ACTN-T-TNNC", "")
	exp.AddSequence("C", "ACTN-TNT--N", "")
	exp.AddSequence("D", "N-ANNNNNNCC", "")

	exp2 = NewAlign(NUCLEOTIDS)
	exp2.AddSequence("A", "ACANGA-TACC", "")
	exp2.AddSequence("B", "ACTN-T-TNNC", "")
	exp2.AddSequence("C", "ACTN-TNT--N", "")
	exp2.AddSequence("D", "N-ANNNNNNCC", "")

	exp3 = NewAlign(NUCLEOTIDS)
	exp3.AddSequence("A", "ACANGA-TACC", "")
	exp3.AddSequence("B", "ACNN-N-TNNC", "")
	exp3.AddSequence("C", "ACNN-NNT--N", "")
	exp3.AddSequence("D", "N-ANNNNNNCC", "")

	res, _ = a.Clone()
	res.MaskUnique("", "AMBIG")

	if !res.Identical(exp) {
		t.Error(fmt.Errorf("Mask unique failed"))
	}

	res, _ = a.Clone()
	res.MaskUnique("A", "AMBIG")
	if !res.Identical(exp2) {
		t.Error(fmt.Errorf("Mask unique with refseq failed"))
	}

	res, _ = a.Clone()
	res.MaskOccurences("A", 2, "AMBIG")
	if !res.Identical(exp3) {
		t.Error(fmt.Errorf("Mask unique with refseq and maxocc failed"))
	}
}

func Test_align_MaskUniqueGAP(t *testing.T) {
	var a, res, exp, exp2, exp3 Alignment

	a = NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "ACANGA-TACC", "")
	a.AddSequence("B", "ACTN-T-TTTC", "")
	a.AddSequence("C", "ACTN-TTT--T", "")
	a.AddSequence("D", "C-ANCCCCCCC", "")

	exp = NewAlign(NUCLEOTIDS)
	exp.AddSequence("A", "ACAN---T-CC", "")
	exp.AddSequence("B", "ACTN-T-T--C", "")
	exp.AddSequence("C", "ACTN-T-T---", "")
	exp.AddSequence("D", "--AN-----CC", "")

	exp2 = NewAlign(NUCLEOTIDS)
	exp2.AddSequence("A", "ACANGA-TACC", "")
	exp2.AddSequence("B", "ACTN-T-T--C", "")
	exp2.AddSequence("C", "ACTN-T-T---", "")
	exp2.AddSequence("D", "--AN-----CC", "")

	exp3 = NewAlign(NUCLEOTIDS)
	exp3.AddSequence("A", "ACANGA-TACC", "")
	exp3.AddSequence("B", "AC-N---T--C", "")
	exp3.AddSequence("C", "AC-N---T---", "")
	exp3.AddSequence("D", "--AN-----CC", "")

	res, _ = a.Clone()
	res.MaskUnique("", "GAP")

	if !res.Identical(exp) {
		t.Error(fmt.Errorf("Mask unique failed"))
	}

	res, _ = a.Clone()
	res.MaskUnique("A", "GAP")
	if !res.Identical(exp2) {
		t.Error(fmt.Errorf("Mask unique with refseq failed"))
	}

	res, _ = a.Clone()
	res.MaskOccurences("A", 2, "GAP")
	if !res.Identical(exp3) {
		t.Error(fmt.Errorf("Mask unique with refseq and maxocc failed"))
	}
}

func Test_align_MaskUniqueMAJ(t *testing.T) {
	var a, res, exp, exp2, exp3 Alignment

	a = NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "AAAAAAAAAAA", "")
	a.AddSequence("B", "AAAAAAAAAAA", "")
	a.AddSequence("C", "CCCCCCCCCCC", "")
	a.AddSequence("D", "CCCCCCCCCCC", "")
	a.AddSequence("E", "TTTTTTTTTTT", "")

	exp = NewAlign(NUCLEOTIDS)
	exp.AddSequence("A", "AAAAAAAAAAA", "")
	exp.AddSequence("B", "AAAAAAAAAAA", "")
	exp.AddSequence("C", "CCCCCCCCCCC", "")
	exp.AddSequence("D", "CCCCCCCCCCC", "")
	exp.AddSequence("E", "AAAAAAAAAAA", "")

	exp2 = NewAlign(NUCLEOTIDS)
	exp2.AddSequence("A", "AAAAAAAAAAA", "")
	exp2.AddSequence("B", "AAAAAAAAAAA", "")
	exp2.AddSequence("C", "CCCCCCCCCCC", "")
	exp2.AddSequence("D", "CCCCCCCCCCC", "")
	exp2.AddSequence("E", "CCCCCCCCCCC", "")

	exp3 = NewAlign(NUCLEOTIDS)
	exp3.AddSequence("A", "AAAAAAAAAAA", "")
	exp3.AddSequence("B", "AAAAAAAAAAA", "")
	exp3.AddSequence("C", "CCCCCCCCCCC", "")
	exp3.AddSequence("D", "CCCCCCCCCCC", "")
	exp3.AddSequence("E", "CCCCCCCCCCC", "")

	res, _ = a.Clone()
	res.MaskUnique("", "MAJ")

	if !res.Identical(exp) {
		t.Error(fmt.Errorf("Mask unique failed"))
	}

	res, _ = a.Clone()
	res.MaskUnique("A", "MAJ")
	if !res.Identical(exp2) {
		fmt.Println(res)
		fmt.Println(exp2)
		t.Error(fmt.Errorf("Mask unique with refseq failed"))
	}

	res, _ = a.Clone()
	res.MaskOccurences("A", 2, "MAJ")
	if !res.Identical(exp3) {
		t.Error(fmt.Errorf("Mask unique with refseq and maxocc failed"))
	}
}

func Test_align_InverseCoordinates(t *testing.T) {
	var a Alignment = NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "ACANGA-TACC", "")
	a.AddSequence("B", "ACTN-T-TTTC", "")
	a.AddSequence("C", "ACTN-TTT--T", "")
	a.AddSequence("D", "C-ANCCCCCCC", "")

	expinvlengths := []int{2, 3}
	expinvstarts := []int{0, 8}

	invstarts, invlengths, err := a.InverseCoordinates(0, 11)
	if err != nil {
		t.Error(err)
	}
	if len(invstarts) != 0 {
		t.Errorf("Invstarts should be 0 length")
	}
	if len(invlengths) != 0 {
		t.Errorf("Invlengths should be 0 length")
	}

	invstarts, invlengths, err = a.InverseCoordinates(2, 6)
	if err != nil {
		t.Error(err)
	}

	if len(invstarts) != 2 {
		t.Errorf("Invstarts should be 2 length")
	}
	if len(invlengths) != 2 {
		t.Errorf("Invlengths should be 2 length")
	}

	if !reflect.DeepEqual(expinvlengths, invlengths) {
		t.Error(fmt.Errorf("Invlengths is not expected, have %v, want %v", invlengths, expinvlengths))
	}
	if !reflect.DeepEqual(expinvstarts, invstarts) {
		t.Error(fmt.Errorf("Invlengths is not expected, have %v, want %v", invlengths, expinvlengths))
	}
}

func Test_align_InversePositions(t *testing.T) {
	var a Alignment = NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "ACANGA-TACC", "")
	a.AddSequence("B", "ACTN-T-TTTC", "")
	a.AddSequence("C", "ACTN-TTT--T", "")
	a.AddSequence("D", "C-ANCCCCCCC", "")

	positions := []int{1, 2, 6, 10}
	expinvpositions := []int{0, 3, 4, 5, 7, 8, 9}

	invpositions, err := a.InversePositions(positions)
	if err != nil {
		t.Error(err)
	}
	if !reflect.DeepEqual(expinvpositions, invpositions) {
		t.Error(fmt.Errorf("Invpositions is not expected, have %v, want %v", invpositions, expinvpositions))
	}

	positions = []int{}
	expinvpositions = []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}

	invpositions, err = a.InversePositions(positions)
	if err != nil {
		t.Error(err)
	}
	if !reflect.DeepEqual(expinvpositions, invpositions) {
		t.Error(fmt.Errorf("Invpositions is not expected, have %v, want %v", invpositions, expinvpositions))
	}

	positions = []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}
	expinvpositions = []int{}

	invpositions, err = a.InversePositions(positions)
	if err != nil {
		t.Error(err)
	}
	if !reflect.DeepEqual(expinvpositions, invpositions) {
		t.Error(fmt.Errorf("Invpositions is not expected, have %v, want %v", invpositions, expinvpositions))
	}
}

func Test_align_RefSites(t *testing.T) {
	var err error
	var alisites []int

	in := NewAlign(UNKNOWN)

	in.AddSequence("Seq0000", "--ACG--AT---GC", "")
	in.AddSequence("Seq0001", "GGACGTTATCGGGC", "")
	in.AutoAlphabet()

	expalisites := []int{2, 3, 4, 7}
	sites := []int{0, 1, 2, 3}

	alisites, err = in.RefSites("Seq0000", sites)
	if err != nil {
		t.Error(err)
	}
	if !reflect.DeepEqual(expalisites, alisites) {
		t.Error(fmt.Errorf("alisites: expected: %v, have: %v", expalisites, alisites))
	}
}

func Test_align_InformativeSites(t *testing.T) {
	atab := make([]*align, 0)
	a := NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "AAAAAAAAAAA", "")
	a.AddSequence("B", "AAAAAAAAAAA", "")
	a.AddSequence("C", "AAAAAAAAAAA", "")
	a.AddSequence("D", "AAAAAAAAAAA", "")
	atab = append(atab, a)

	a = NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "ACGTACGTACG", "")
	atab = append(atab, a)

	a = NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "AAAAAAAAAAA", "")
	a.AddSequence("B", "AAAAAAAAAAA", "")
	a.AddSequence("C", "AAAACCAAAAA", "")
	a.AddSequence("D", "CCCCCCCCCCC", "")
	atab = append(atab, a)

	a = NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "AAA-A--AAAA", "")
	a.AddSequence("B", "AAAAAAA--AA", "")
	a.AddSequence("C", "AAAACCAAA--", "")
	a.AddSequence("D", "CCC--CCCCCC", "")
	atab = append(atab, a)

	a = NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "AAAAA--AAAA", "")
	a.AddSequence("B", "AAAAAAA--AA", "")
	a.AddSequence("C", "AAAACCAAA--", "")
	a.AddSequence("D", "CCC--CCCCCC", "")
	a.AddSequence("E", "CCC-CCCCCCC", "")
	atab = append(atab, a)

	a = NewAlign(NUCLEOTIDS)
	a.AddSequence("A", "AAAAA--AAAA", "")
	a.AddSequence("B", "AAAAAAA--AA", "")
	a.AddSequence("C", "AAAACCAAA--", "")
	a.AddSequence("D", "CCC--CCCCCC", "")
	a.AddSequence("E", "CCC-CCCCCCC", "")
	a.AddSequence("D", "GCC--CCCCCC", "")
	a.AddSequence("E", "GCC-CCCCCCC", "")
	atab = append(atab, a)

	tests := []struct {
		name      string
		fields    *align
		wantSites []int
	}{
		{name: "Identical", fields: atab[0], wantSites: []int{}},
		{name: "Single sequence", fields: atab[1], wantSites: []int{}},
		{name: "2 informative", fields: atab[2], wantSites: []int{4, 5}},
		{name: "Gaps no informative", fields: atab[3], wantSites: []int{}},
		{name: "Gaps 9 informatives", fields: atab[4], wantSites: []int{0, 1, 2, 4, 6, 7, 8, 9, 10}},
		{name: "Gaps 9 informatives /2", fields: atab[4], wantSites: []int{0, 1, 2, 4, 6, 7, 8, 9, 10}},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotSites := tt.fields.InformativeSites(); !reflect.DeepEqual(gotSites, tt.wantSites) {
				t.Errorf("align.InformativeSites() = %v, want %v", gotSites, tt.wantSites)
			}
		})
	}
}
