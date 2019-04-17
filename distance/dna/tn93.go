package dna

import (
	"math"

	"github.com/evolbioinfo/goalign/align"
)

type TN93Model struct {
	/* Vector of nt proba */
	pi            []float64 // proba of each nt
	numSites      float64   // Number of selected sites (no gaps)
	selectedSites []bool    // true for selected sites
	removegaps    bool      // If true, we will remove posision with >=1 gaps
	gamma         bool
	alpha         float64
}

func NewTN93Model(removegaps bool) *TN93Model {
	return &TN93Model{
		nil,
		0,
		nil,
		removegaps,
		false,
		0.,
	}
}

/* computes TN93 distance between 2 sequences */
func (m *TN93Model) Distance(seq1 []rune, seq2 []rune, weights []float64) (float64, error) {
	var dist float64
	var b1, b2, b3 float64

	trS, trV, p1, p2, total := countMutations(seq1, seq2, m.selectedSites, weights)
	trS, trV, p1, p2 = trS/total, trV/total, p1/total, p2/total

	piy := m.pi[1] + m.pi[3]
	pir := m.pi[0] + m.pi[2]
	papg := m.pi[0] * m.pi[2]
	pcpt := m.pi[1] * m.pi[3]

	y := papg / (papg + pcpt)
	e1 := 1 - trV/(2*piy*pir)
	e2 := 1 - trV/(2*pir) - pir*p1/(2*papg)
	e3 := 1 - trV/(2*piy) - piy*p2/(2*pcpt)

	if m.gamma {
		b1 = (piy/pir*m.alpha*(1.-math.Pow(e1, -1./m.alpha)) - 1./pir*m.alpha*(1.-math.Pow(e2, -1./m.alpha)))
		b2 = (pir/piy*m.alpha*(1.-math.Pow(e1, -1./m.alpha)) - 1./piy*m.alpha*(1.-math.Pow(e3, -1./m.alpha)))
		b3 = -m.alpha * (1. - math.Pow(e1, -1./m.alpha))
	} else {
		b1 = piy/pir*math.Log(e1) - 1./pir*math.Log(e2)
		b2 = pir/piy*math.Log(e1) - 1./piy*math.Log(e3)
		b3 = -math.Log(e1)
	}

	dist = 2.*(m.pi[0]*m.pi[2]+m.pi[1]*m.pi[3])*(y*b1+(1-y)*b2) + 2*pir*piy*b3
	if dist > 0 {
		return dist, nil
	} else {
		return 0, nil
	}
}

func (m *TN93Model) InitModel(al align.Alignment, weights []float64, gamma bool, alpha float64) (err error) {
	m.gamma = gamma
	m.alpha = alpha
	m.numSites, m.selectedSites = selectedSites(al, weights, m.removegaps)
	m.pi, err = probaNt(al, m.selectedSites, weights)
	return
}
