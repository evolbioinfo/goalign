package distance

import (
	"fmt"
	"github.com/fredericlemoine/goalign/align"
	"math"
)

type F81Model struct {
	/* Vector of nt proba */
	pi            []float64 // proba of each nt
	b1            float64   // 1- sum(pi^2)
	numSites      float64   // Number of selected sites (no gaps)
	selectedSites []bool    // true for selected sites
}

func NewF81Model() *F81Model {
	return &F81Model{
		make([]float64, 4),
		0,
		0,
		nil,
	}
}

/* computes F81 distance between 2 sequences */
func (m *F81Model) Distance(seq1 []rune, seq2 []rune, weights []float64) float64 {
	diff, total := countDiffs(seq1, seq2, weights)
	diff = diff / total
	dist := -1.0 * m.b1 * math.Log(1.0-diff/m.b1)
	if dist > 0 {
		return (dist)
	} else {
		return (0)
	}
}

func (m *F81Model) InitModel(al align.Alignment, weights []float64) {
	m.numSites, m.selectedSites = selectedSites(al, weights)
	total := 0.0
	for i := 0; i < al.NbSequences(); i++ {
		seq1, _ := al.GetSequenceChar(i)
		for _, n := range seq1 {
			if m.selectedSites[n] {
				switch n {
				case 'A':
					m.pi[0]++
				case 'C':
					m.pi[1]++
				case 'G':
					m.pi[2]++
				case 'T':
					m.pi[3]++
				}
				total++
			}
		}
	}
	for i, _ := range m.pi {
		m.pi[i] /= total
		m.b1 += math.Pow(m.pi[i], 2.0)
	}
	m.b1 = 1 - m.b1
	fmt.Println(m.b1)
}
