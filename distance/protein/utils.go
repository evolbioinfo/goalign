package protein

import (
	"fmt"
	"math"

	"github.com/evolbioinfo/goalign/align"
)

type seqpairdist struct {
	i, j       int
	seq1, seq2 []rune
	seq1Ambigu []bool
	seq2Ambigu []bool
}

/* Returns the sites of the alignments that contains only nucleotides and no gaps */
func selectedSites(al align.Alignment, weights []float64, removeGappedPositions bool) (float64, []bool) {
	selectedSites := make([]bool, al.Length())
	numSites := 0.0
	for l := 0; l < al.Length(); l++ {
		w := 1.0
		if weights != nil {
			w = weights[l]
		}
		selectedSites[l] = true
		for i := 0; i < al.NbSequences() && removeGappedPositions; i++ {
			seq, _ := al.GetSequenceCharById(i)
			if al.AlphabetCharToIndex(seq[l]) == -1 || seq[l] == '*' || seq[l] == '?' || seq[l] == '-' {
				selectedSites[l] = false
			}
		}
		if selectedSites[l] {
			numSites += w
		}
	}
	return numSites, selectedSites
}

func aaFrequency(a align.Alignment, weights []float64, selected []bool) ([]float64, error) {
	var i, j int
	var w, sum float64
	if a.Alphabet() != align.AMINOACIDS {
		return nil, fmt.Errorf("Alphabet is not AminoAcids")
	}
	ns := len(a.AlphabetCharacters())
	var num []float64 = make([]float64, 20)
	var freq []float64 = make([]float64, 20)

	for i = 0; i < ns; i++ {
		freq[i] = 1. / float64(ns)
		num[i] = 0.0
	}

	// Count occurences of different amino acids
	a.IterateChar(func(name string, sequence []rune) bool {
		for j = 0; j < len(sequence); j++ {
			if selected[j] {
				w = weights[j]
				idx := a.AlphabetCharToIndex(sequence[j])
				if idx >= 0 {
					num[idx] += w
				} else {
					for i = 0; i < ns; i++ {
						num[i] = w * freq[i]
					}
				}
			}
		}
		return false
	})

	// if at least one frequency equals 0 then add a pseudo-count
	// as these are doubles, cannot test equality to 0, then test less than minimum value it can have (1./20)
	oneLessThanCutoff := false
	for _, v := range num {
		if v < 1./float64(ns) {
			oneLessThanCutoff = true
			break
		}
	}

	for i, v := range num {
		if oneLessThanCutoff {
			num[i] = v + 1.0
		}
		sum += num[i]
	}
	for i, _ := range num {
		freq[i] = num[i] / sum
	}

	return freq, nil
}

func isAmbigu(c rune) bool {
	return (c == align.GAP || c == align.POINT || c == align.OTHER || c == align.ALL_AMINO)
}

func sign(a, b float64) float64 {
	if b > .0 {
		return math.Abs(a)
	} else {
		return -math.Abs(a)
	}
}

func shift(a, b, c, d *float64) {
	(*a) = (*b)
	(*b) = (*c)
	(*c) = (*d)
}

func checkAmbiguities(pair *seqpairdist, stepsize int) (ambig []bool) {
	var j int
	pair.seq1Ambigu = make([]bool, len(pair.seq1))
	pair.seq2Ambigu = make([]bool, len(pair.seq2))

	for j = 0; j < len(pair.seq1); j += stepsize {
		pair.seq1Ambigu[j] = false
		pair.seq2Ambigu[j] = false
		if isAmbigu(pair.seq1[j]) {
			pair.seq1Ambigu[j] = true
		}
		if isAmbigu(pair.seq2[j]) {
			pair.seq2Ambigu[j] = true
		}
	}
	return
}

func check2SequencesDiff(pair *seqpairdist) (ret bool) {
	var i int
	for i = 0; i < len(pair.seq1); i++ {
		if (!pair.seq1Ambigu[i] && !pair.seq2Ambigu[i]) && (pair.seq1[i] != pair.seq2[i]) {
			return true
		}
	}
	return false
}
