package dna

import (
	"math"

	"github.com/evolbioinfo/goalign/align"
	"gonum.org/v1/gonum/mat"
)

type TN93Model struct {
	/* Vector of nt proba */
	pi            []float64 // proba of each nt
	numSites      float64   // Number of selected sites (no gaps)
	selectedSites []bool    // true for selected sites
	removegaps    bool      // If true, we will remove posision with >=1 gaps
	// Parameters (for eigen values/vectors computation)
	// See https://en.wikipedia.org/wiki/Models_of_DNA_evolution#F81_model_(Felsenstein_1981)
	qmatrix mat.Dense
}

func NewTN93Model(removegaps bool) *TN93Model {
	return &TN93Model{
		nil,
		0,
		nil,
		removegaps,
	}
}

/* computes TN93 distance between 2 sequences */
func (m *TN93Model) Distance(seq1 []rune, seq2 []rune, weights []float64) (float64, error) {
	trS, trV, p1, p2, total := countMutations(seq1, seq2, m.selectedSites, weights)
	trS, trV, p1, p2 = trS/total, trV/total, p1/total, p2/total

	piy := m.pi[1] + m.pi[3]
	pir := m.pi[0] + m.pi[2]
	y := m.pi[0] * m.pi[2] / (m.pi[0]*m.pi[2] + m.pi[1]*m.pi[3])
	e1 := 1 - trV/(2*piy*pir)
	e2 := 1 - trV/(2*pir) - pir*p1/(2*m.pi[0]*m.pi[2])
	e3 := 1 - trV/(2*piy) - piy*p2/(2*m.pi[1]*m.pi[3])
	b1 := piy/pir*math.Log(e1) - 1/pir*math.Log(e2)
	b2 := pir/piy*math.Log(e1) - 1/piy*math.Log(e3)
	b3 := -math.Log(e1)

	dist := 2*(m.pi[0]*m.pi[2]+m.pi[1]*m.pi[3])*(y*b1+(1-y)*b2) + 2*pir*piy*b3
	if dist > 0 {
		return dist, nil
	} else {
		return 0, nil
	}
}

func (m *TN93Model) InitModel(al align.Alignment, weights []float64) (err error) {
	m.numSites, m.selectedSites = selectedSites(al, weights, m.removegaps)
	m.pi, err = probaNt(al, m.selectedSites, weights)
	return
}

func (m *F81Model) SetParameters(kappa1, kappa2, piA, piC, piG, piT float64) {
	m.qmatrix = mat.NewDense(4, 4, []float64{
		-(piC + kappa1*piG + piT), piC, kappa1 * piG, piT,
		piA, -(piA + piG + kappa2*piT), piG, kappa2 * piT,
		kappa1 * piA, piC, -(kappa1*piA + piC + piT), piT,
		piA, kappa2 * piC, piG, -(piA + kappa2*piC + piG),
	})
	// Normalization of Q
	norm := 1. / (2. * (piA*piC + piC*piG + piA*piT + piG*piT + kappa2*piC*piT + kappa1*piA*piG))
	m.qmatrix.Apply(func(i, j int, v float64) float64 { v * norm }, m.qmatrix)
}

func (m *F81Model) Eigens() (val []float64, leftvector, rightvector [][]float64, err error) {
	// Compute eigen values, left and right eigenvectors of Q
	eigen := &mat.Eigen{}
	if ok = eigen.Factorize(model.qmatrix, mat.EigenRight); !ok {
		err = fmt.Errorf("Problem during matrix decomposition")
		return
	}

	val = make([]float64, 4)
	for i, b := range model.eigen.Values(nil) {
		val[i] = real(b)
	}
	u := model.eigen.VectorsTo(nil)
	reigenvect := mat.NewDense(4, 4, nil)
	leigenvect := mat.NewDense(4, 4, nil)
	reigenvect.Apply(func(i, j int, val float64) float64 { return real(u.At(i, j)) }, reigenvect)
	leigenvect.Inverse(model.reigenvect)

	leftvector = [][]float64{
		leigenvect.RawRowView(0),
		leigenvect.RawRowView(1),
		leigenvect.RawRowView(2),
		leigenvect.RawRowView(3),
	}
	rightvector = [][]float64{
		reigenvect.RawRowView(0),
		reigenvect.RawRowView(1),
		reigenvect.RawRowView(2),
		reigenvect.RawRowView(3),
	}

	return
}
