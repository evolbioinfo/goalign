package dna

import (
	"gonum.org/v1/gonum/mat"
)

const (
	DBL_MIN = 2.2250738585072014e-308
)

type DNAModel interface {
	Eigens() (val []float64, leftvectors, rightvectors *mat.Dense, err error)
	Analytical() bool                // returns true if analytical pij computation is possible and implemented
	Pij(i, j int, l float64) float64 // Returns -1 if not possible to compute it anatically without eigens (or not yet implemented)
}
