package dna

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"gonum.org/v1/gonum/mat"
)

type GTRModel struct {
	qmatrix *mat.Dense
}

func NewGTRModel() *GTRModel {
	return &GTRModel{
		nil,
	}
}

/* computes TN82 distance between 2 sequences */
func (m *GTRModel) Distance(seq1 []rune, seq2 []rune, weights []float64) (dist float64, err error) {
	err = fmt.Errorf("Cannot compute distances on GTR model")
	return
}

func (m *GTRModel) InitModel(al align.Alignment, weights []float64) (err error) {

	return
}

//  /          \
// | *  d  -  b |
// | d  *  e  a |
// | -  e  *  c |
// | b  a  c  * |
//  \          /
func (m *GTRModel) SetParameters(a, b, c, d, e, piA, piC, piG, piT float64) {
	m.qmatrix = mat.NewDense(4, 4, []float64{
		-(d*piC + piG + b*piT), d * piC, piG, b * piT,
		d * piA, -(d*piA + e*piG + a*piT), e * piG, a * piT,
		piA, e * piC, -(piA + e*piC + c*piT), c * piT,
		b * piA, a * piC, c * piG, -(b*piA + a*piC + c*piG),
	})
	// Normalization of Q
	norm := 1. / (piA*(d+1+b) + piC*(d+e+a) + piG*(1+e+c) + piT*(b+a+c))
	m.qmatrix.Apply(func(i, j int, v float64) float64 { return v * norm }, m.qmatrix)
}

func (m *GTRModel) Eigens() (val []float64, leftvector, rightvector [][]float64, err error) {
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
