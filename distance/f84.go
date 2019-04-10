package distance

import (
	"math"

	"github.com/evolbioinfo/goalign/align"
)

type F84Model struct {
	pi            []float64 // Vector of nt stationary proba
	a, b, c       float64   // Parameters for distance computation
	numSites      float64   // Number of selected sites (no gaps)
	selectedSites []bool    // true for selected sites
	removegaps    bool      // If true, we will remove posision with >=1 gaps
	// Parameters (for eigen values/vectors computation)
	// https://en.wikipedia.org/wiki/Models_of_DNA_evolution#HKY85_model_(Hasegawa,_Kishino_and_Yano_1985)
	piA, piC, piG, piT float64
	kappa              float64
}

func NewF84Model(removegaps bool) *F84Model {
	return &F84Model{
		nil,
		0, 0, 0,
		0,
		nil,
		removegaps,
		1.0,
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

func (m *F81Model) SetParameters(kappa, piA, piC, piG, piT float64) {
	//m.qmatrix = mat.NewDense(4, 4, []float64{
	//	-(piC + (1+kappa/piR)*piG + piT), piC, (1 + kappa/piR) * piG, piT,
	//	piA, -(piA + piG + (1+kappa/piY)*piT), piG, (1 + kappa/piY) * piT,
	//	(1 + kappa/piR) * piA, piC, -((1+kappa/piR)*piA + piC + piT), piT,
	//	piA, (1 + kappa/piY) * piC, piG, -(piA + (1+kappa/piY)*piC + piG),
	//})
	// Normalization of Q
	m.kappa = kappa
	m.piA = piA
	m.piC = piC
	m.piG = piG
	m.piT = piT
}

// See http://biopp.univ-montp2.fr/Documents/ClassDocumentation/bpp-phyl/html/F84_8cpp_source.html
func (m *F81Model) Eigens() (val []float64, leftvector, rightvector [][]float64, err error) {
	piY := m.piT + m.piC
	piR := m.piA + m.piG
	norm := 1. / (1 - m.piA_*m.piA_ - m.piC_*m.piC_ - m.piG_*m.piG_ - m.piT_*m.piT_ + 2.*m.kappa_*(m.piC_*m.piT_/piY_+m.piA_*m.piG_/piR_))

	val = []float64{
		0,
		-norm * (1 + m.kappa),
		-norm * (1 + m.kappa),
		-norm,
	}

	leftvector = [][]float64{
		[]float64{m.piA, m.piC, m.piG, m.piT},
		[]float64{0., m.piT / piY, 0., -m.piT / piY},
		[]float64{m.piG / piR, 0., -piG / piR, 0.},
		[]float64{m.piA * piY / piR, -piC, m.piG * piY / piR, -m.piT},
	}

	rightvector = [][]float64{
		[]float64{1., 0., 1., 1.},
		[]float64{1., 1., 0., -piR / piY},
		[]float64{1., 0., -m.piA / m.piG, 1.},
		[]float64{1., -m.piC / m.piT, 0., -piR / piY},
	}

	return
}
