package dna

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

const (
	DBL_MIN = 2.2250738585072014e-308
)

type DNAModel interface {
	Eigens() (val []float64, leftvectors, rightvectors []float64, err error)
}

// Probability matrix
type Pij struct {
	length float64  // branch length / t
	model  DNAModel // Model
	pij    *mat.Dense
}

func NewPij(m DNAModel, l float64) (pij *Pij, err error) {
	pij = &Pij{l, m, mat.NewDense(4, 4, nil)}
	err = pij.SetLength(l)
	return
}

func (pij *Pij) SetLength(l float64) (err error) {
	var i int
	var v []float64
	var left, right []float64
	var expt []float64
	var uexpt *mat.Dense
	ns := 4

	if v, left, right, err = pij.model.Eigens(); err != nil {
		return
	}

	expt = make([]float64, ns)
	for i = 0; i < ns; i++ {
		expt[i] = float64(math.Exp(v[i] * l))
	}

	uexpt = mat.NewDense(ns, ns, nil)
	lMat := mat.NewDense(ns, ns, left).T()
	rMat := mat.NewDense(ns, ns, right).T()

	uexpt.Mul(lMat, mat.NewDiagDense(ns, expt))
	pij.pij.Mul(uexpt, rMat)

	return
}

func (pij *Pij) Pij(i, j int) float64 {
	return pij.pij.At(i, j)
}
