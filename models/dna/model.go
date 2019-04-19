package dna

import (
	"math"

	"gonum.org/v1/gonum/mat"
)

const (
	DBL_MIN = 2.2250738585072014e-308
)

type DNAModel interface {
	Eigens() (val []float64, leftvectors, rightvectors *mat.Dense, err error)
	analytical() bool                // returns true if analytical pij computation is possible and implemented
	pij(i, j int, l float64) float64 // Returns -1 if not possible to compute it anatically without eigens (or not yet implemented)
}

// Probability matrix
type Pij struct {
	length float64  // branch length / t
	model  DNAModel // Model
	pij    *mat.Dense
	expt   []float64  // tmp array
	uexpt  *mat.Dense // tmp Dense
}

func NewPij(m DNAModel, l float64) (pij *Pij, err error) {
	pij = &Pij{DBL_MIN, m, mat.NewDense(4, 4, []float64{1.0, 0., 0., 0., 0., 1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 1.}), make([]float64, 4), mat.NewDense(4, 4, nil)}
	err = pij.SetLength(l)
	return
}

func (pij *Pij) SetLength(l float64) (err error) {
	if pij.length != l && !pij.model.analytical() {
		var i, j, k int
		var val []float64
		var left, right *mat.Dense
		ns := 4

		if val, left, right, err = pij.model.Eigens(); err != nil {
			return
		}

		for i = 0; i < ns; i++ {
			pij.expt[i] = float64(math.Exp(val[i] * l))
		}
		for i = 0; i < ns; i++ {
			for j = 0; j < ns; j++ {
				pij.uexpt.Set(i, j, right.At(i, j)*pij.expt[j])
			}
		}
		v := 0.0
		for i = 0; i < ns; i++ {
			for j = 0; j < ns; j++ {
				v = 0.0
				for k = 0; k < ns; k++ {
					v += pij.uexpt.At(i, k) * left.At(k, j)
				}
				if v < DBL_MIN {
					v = DBL_MIN
				}
				pij.pij.Set(i, j, v)
			}
		}
	}
	pij.length = l

	return
}

func (pij *Pij) Pij(i, j int) float64 {
	if pij.model.analytical() {
		return (pij.model.pij(i, j, pij.length))
	}
	return pij.pij.At(i, j)
}
