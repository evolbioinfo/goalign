package dna

import (
	"fmt"
	"math"

	"github.com/evolbioinfo/goalign/align"
	"gonum.org/v1/gonum/mat"
)

type F81Model struct {
	pi            []float64 // Vector of nt stationary proba
	b1            float64   // Parameter for distance computation
	numSites      float64   // Number of selected sites (no gaps)
	selectedSites []bool    // true for selected sites
	removegaps    bool      // If true, we will remove posision with >=1 gaps

	// Parameters (for eigen values/vectors computation)
	// See https://en.wikipedia.org/wiki/Models_of_DNA_evolution#F81_model_(Felsenstein_1981)
	qmatrix *mat.Dense
}

func NewF81Model(removegaps bool) *F81Model {
	return &F81Model{
		nil,
		0,
		0,
		nil,
		removegaps,
		nil,
	}
}

/* computes F81 distance between 2 sequences */
func (m *F81Model) Distance(seq1 []rune, seq2 []rune, weights []float64) (float64, error) {
	diff, total := countDiffs(seq1, seq2, m.selectedSites, weights)
	diff = diff / total

	dist := -1.0 * m.b1 * math.Log(1.0-diff/m.b1)
	if dist > 0 {
		return dist, nil
	} else {
		return 0, nil
	}
}

func (m *F81Model) InitModel(al align.Alignment, weights []float64) (err error) {
	m.numSites, m.selectedSites = selectedSites(al, weights, m.removegaps)
	m.b1 = 0.0
	m.pi, err = probaNt(al, m.selectedSites, weights)
	if err == nil {
		for i, _ := range m.pi {
			m.b1 += m.pi[i] * m.pi[i]
		}
		m.b1 = 1 - m.b1
	}
	return
}

func (m *F81Model) SetParameters(piA, piC, piG, piT float64) {
	m.qmatrix = mat.NewDense(4, 4, []float64{
		-(piC + piG + piT), piC, piG, piT,
		piA, -(piA + piG + piT), piG, piT,
		piA, piC, -(piA + piC + piT), piT,
		piA, piC, piG, -(piA + piC + piG),
	})
	// Normalization of Q
	norm := 1. / (2 * (piA*piC + piA*piG + piA*piT + piC*piG + piC*piT + piG*piT))
	m.qmatrix.Apply(func(i, j int, v float64) float64 { return v * norm }, m.qmatrix)
}

func (m *F81Model) Eigens() (val []float64, leftvectors, rightvectors []float64, err error) {
	// Compute eigen values, left and right eigenvectors of Q
	eigen := &mat.Eigen{}
	if ok := eigen.Factorize(m.qmatrix, mat.EigenRight); !ok {
		err = fmt.Errorf("Problem during matrix decomposition")
		return
	}

	val = make([]float64, 4)
	for i, b := range eigen.Values(nil) {
		val[i] = real(b)
	}
	u := eigen.VectorsTo(nil)
	reigenvect := mat.NewDense(4, 4, nil)
	leigenvect := mat.NewDense(4, 4, nil)
	reigenvect.Apply(func(i, j int, val float64) float64 { return real(u.At(i, j)) }, reigenvect)
	leigenvect.Inverse(reigenvect)

	leftvectors = leigenvect.RawMatrix().Data
	rightvectors = reigenvect.RawMatrix().Data

	return
}
