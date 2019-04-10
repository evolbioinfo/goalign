package distance

import (
	"github.com/evolbioinfo/goalign/align"
	"math"
)

type F84Model struct {
	pi            []float64 // Vector of nt stationary proba
	a, b, c       float64   // Parameters for distance computation
	numSites      float64   // Number of selected sites (no gaps)
	selectedSites []bool    // true for selected sites
	removegaps    bool      // If true, we will remove posision with >=1 gaps
}

func NewF84Model(removegaps bool) *F84Model {
	return &F84Model{
		nil,
		0, 0, 0,
		0,
		nil,
		removegaps,
	}
}

/* computes F84 distance between 2 sequences */
func (m *F84Model) Distance(seq1 []rune, seq2 []rune, weights []float64) (float64, error) {
	trS, trV, _, _, total := countMutations(seq1, seq2, m.selectedSites, weights)
	trS, trV = trS/total, trV/total
	dist := -2.0*m.a*math.Log(1.0-trS/(2.0*m.a)-(m.a-m.b)*trV/(2.0*m.a*m.c)) + 2.0*(m.a-m.b-m.c)*math.Log(1-trV/(2.0*m.c))
	if dist > 0 {
		return dist, nil
	} else {
		return 0, nil
	}
}

func (m *F84Model) InitModel(al align.Alignment, weights []float64) (err error) {
	m.numSites, m.selectedSites = selectedSites(al, weights, m.removegaps)
	m.pi, err = probaNt(al, m.selectedSites, weights)
	if err == nil {
		m.a = m.pi[0]*m.pi[2]/(m.pi[0]+m.pi[2]) + m.pi[1]*m.pi[3]/(m.pi[1]+m.pi[3])
		m.b = m.pi[0]*m.pi[2] + m.pi[1]*m.pi[3]
		m.c = (m.pi[0] + m.pi[2]) * (m.pi[1] + m.pi[3])
	}
	return
}
