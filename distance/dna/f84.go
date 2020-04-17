package dna

import (
	"fmt"
	"math"

	"github.com/evolbioinfo/goalign/align"
)

type F84Model struct {
	pi            []float64 // Vector of nt stationary proba
	a, b, c       float64   // Parameters for distance computation
	numSites      float64   // Number of selected sites (no gaps)
	selectedSites []bool    // true for selected sites
	removegaps    bool      // If true, we will remove posision with >=1 gaps
	gamma         bool
	alpha         float64
	sequenceCodes [][]int // Sequences converted to codes
}

func NewF84Model(removegaps bool) *F84Model {
	return &F84Model{
		nil,
		0, 0, 0,
		0,
		nil,
		removegaps,
		false,
		0.,
		nil,
	}
}

/* computes F84 distance between 2 sequences */
func (m *F84Model) Distance(seq1 []int, seq2 []int, weights []float64) (float64, error) {
	var dist float64

	trS, trV, _, _, total := countMutations(seq1, seq2, m.selectedSites, weights)
	trS, trV = trS/total, trV/total
	if m.gamma {
		dist = 2.0 * m.alpha * (m.a*math.Pow((1.0-trS/(2.0*m.a)-(m.a-m.b)*trV/(2.0*m.a*m.c)), -1./m.alpha) +
			(m.b+m.c-m.a)*math.Pow((1-trV/(2.0*m.c)), -1./m.alpha) -
			m.b - m.c)
	} else {
		dist = -2.0*m.a*math.Log(1.0-trS/(2.0*m.a)-(m.a-m.b)*trV/(2.0*m.a*m.c)) + 2.0*(m.a-m.b-m.c)*math.Log(1-trV/(2.0*m.c))
	}

	return dist, nil
}

func (m *F84Model) InitModel(al align.Alignment, weights []float64, gamma bool, alpha float64) (err error) {
	m.gamma = gamma
	m.alpha = alpha
	m.numSites, m.selectedSites = selectedSites(al, weights, m.removegaps)
	if m.sequenceCodes, err = alignmentToCodes(al); err != nil {
		return
	}
	m.pi, err = probaNt(m.sequenceCodes, m.selectedSites, weights)
	if err == nil {
		m.a = m.pi[0]*m.pi[2]/(m.pi[0]+m.pi[2]) + m.pi[1]*m.pi[3]/(m.pi[1]+m.pi[3])
		m.b = m.pi[0]*m.pi[2] + m.pi[1]*m.pi[3]
		m.c = (m.pi[0] + m.pi[2]) * (m.pi[1] + m.pi[3])
	}
	return
}

// Sequence returns the ith sequence of the alignment
// encoded in int
func (m *F84Model) Sequence(i int) (seq []int, err error) {
	if i < 0 || i >= len(m.sequenceCodes) {
		err = fmt.Errorf("This sequence does not exist: %d", i)
		return
	}
	seq = m.sequenceCodes[i]
	return
}
