package dna

import (
	"math"
	"testing"

	"github.com/evolbioinfo/goalign/models"
	"gonum.org/v1/gonum/mat"
)

// Tests models with parameters making them equal to JC
const (
	TEST_PRECISION = 0.00001
)

func TestJCEigens(t *testing.T) {
	var v []float64
	var l, r *mat.Dense
	var vMat, lMat, rMat mat.Matrix
	var err error
	var expQMatrix *mat.Dense
	var resQMatrix *mat.Dense = mat.NewDense(4, 4, nil)
	var sub *mat.Dense = mat.NewDense(4, 4, nil)

	expQMatrix = mat.NewDense(4, 4, []float64{
		-1.0000000, 0.3333333, 0.3333333, 0.3333333,
		0.3333333, -1.0000000, 0.3333333, 0.3333333,
		0.3333333, 0.3333333, -1.0000000, 0.3333333,
		0.3333333, 0.3333333, 0.3333333, -1.0000000,
	})

	m := NewJCModel()
	m.InitModel()
	if v, l, r, err = m.Eigens(); err != nil {
		t.Errorf("Error while computing JC eigen vectors: %v", err)
	}

	// Transpose because it was in col-major format
	lMat = l
	rMat = r

	vMat = mat.NewDiagDense(4, v)

	resQMatrix.Mul(rMat, vMat)
	resQMatrix.Mul(resQMatrix, lMat)
	sub.Sub(resQMatrix, expQMatrix)

	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			if math.Abs(sub.At(i, j)) > TEST_PRECISION {
				t.Errorf("Expected QMatrix is different from Resulting QMatrix %v", sub.At(i, j))
			}
		}
	}
}

// func TestJCPij(t *testing.T) {
// 	var v, l, r []float64
// 	var vMat, lMat, rMat mat.Matrix
// 	var err error
// 	var expQMatrix *mat.Dense
// 	var resQMatrix *mat.Dense = mat.NewDense(4, 4, nil)
// 	var sub *mat.Dense = mat.NewDense(4, 4, nil)

// 	m := NewJCModel()
// 	m.InitModel()

// 	for _, l := range []float64{0.1, 0.2, 0.5, 1., 2.} {
// 		pij := NewPij(m, l)
// 		for i := 0; i < 4; i++ {
// 			for j := 0; j < 4; j++ {

// 			}
// 		}
// 	}
// }

func TestF81Eigens(t *testing.T) {
	var v []float64
	var l, r *mat.Dense
	var vMat, lMat, rMat mat.Matrix
	var err error
	var expQMatrix *mat.Dense
	var resQMatrix *mat.Dense = mat.NewDense(4, 4, nil)
	var sub *mat.Dense = mat.NewDense(4, 4, nil)

	expQMatrix = mat.NewDense(4, 4, []float64{
		-1.0000000, 0.3333333, 0.3333333, 0.3333333,
		0.3333333, -1.0000000, 0.3333333, 0.3333333,
		0.3333333, 0.3333333, -1.0000000, 0.3333333,
		0.3333333, 0.3333333, 0.3333333, -1.0000000,
	})
	t.Logf("expQ = %v", mat.Formatted(expQMatrix, mat.Prefix("                 "), mat.Squeeze()))

	m := NewF81Model()
	m.InitModel(1./4., 1./4., 1./4., 1./4.)

	if v, l, r, err = m.Eigens(); err != nil {
		t.Errorf("Error while computing F81 eigen vectors: %v", err)
	}

	// Transpose because it was in col-major formrat
	lMat = l
	rMat = r

	vMat = mat.NewDiagDense(4, v)

	resQMatrix.Mul(rMat, vMat)
	resQMatrix.Mul(resQMatrix, lMat)
	sub.Sub(resQMatrix, expQMatrix)

	t.Logf("resQ = %v", mat.Formatted(resQMatrix, mat.Prefix("                 "), mat.Squeeze()))

	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			if math.Abs(sub.At(i, j)) > TEST_PRECISION {
				t.Errorf("Expected QMatrix is different from Resulting QMatrix %v", sub.At(i, j))
			}
		}
	}
}

func TestF84Eigens(t *testing.T) {
	var v []float64
	var l, r *mat.Dense
	var vMat, lMat, rMat mat.Matrix
	var err error
	var expQMatrix *mat.Dense
	var resQMatrix *mat.Dense = mat.NewDense(4, 4, nil)
	var sub *mat.Dense = mat.NewDense(4, 4, nil)

	expQMatrix = mat.NewDense(4, 4, []float64{
		-1.0000000, 0.3333333, 0.3333333, 0.3333333,
		0.3333333, -1.0000000, 0.3333333, 0.3333333,
		0.3333333, 0.3333333, -1.0000000, 0.3333333,
		0.3333333, 0.3333333, 0.3333333, -1.0000000,
	})
	t.Logf("expQ = %v", mat.Formatted(expQMatrix, mat.Prefix("                 "), mat.Squeeze()))

	m := NewF84Model()
	m.InitModel(0., 1./4., 1./4., 1./4., 1./4.)

	if v, l, r, err = m.Eigens(); err != nil {
		t.Errorf("Error while computing F81 eigen vectors: %v", err)
	}

	// Transpose because it was in col-major formrat
	lMat = l
	rMat = r

	vMat = mat.NewDiagDense(4, v)

	resQMatrix.Mul(rMat, vMat)
	resQMatrix.Mul(resQMatrix, lMat)
	sub.Sub(resQMatrix, expQMatrix)

	t.Logf("resQ = %v", mat.Formatted(resQMatrix, mat.Prefix("                 "), mat.Squeeze()))

	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			if math.Abs(sub.At(i, j)) > TEST_PRECISION {
				t.Errorf("Expected QMatrix is different from Resulting QMatrix %v", sub.At(i, j))
			}
		}
	}
}

func TestGTREigens(t *testing.T) {
	var v []float64
	var l, r *mat.Dense
	var vMat, lMat, rMat mat.Matrix
	var err error
	var expQMatrix *mat.Dense
	var resQMatrix *mat.Dense = mat.NewDense(4, 4, nil)
	var sub *mat.Dense = mat.NewDense(4, 4, nil)

	expQMatrix = mat.NewDense(4, 4, []float64{
		-1.0000000, 0.3333333, 0.3333333, 0.3333333,
		0.3333333, -1.0000000, 0.3333333, 0.3333333,
		0.3333333, 0.3333333, -1.0000000, 0.3333333,
		0.3333333, 0.3333333, 0.3333333, -1.0000000,
	})
	t.Logf("expQ = %v", mat.Formatted(expQMatrix, mat.Prefix("                 "), mat.Squeeze()))

	m := NewGTRModel()
	m.InitModel(1., 1., 1., 1., 1., 1., 1./4., 1./4., 1./4., 1./4.)

	if v, l, r, err = m.Eigens(); err != nil {
		t.Errorf("Error while computing F81 eigen vectors: %v", err)
	}

	// Transpose because it was in col-major formrat
	lMat = l
	rMat = r

	vMat = mat.NewDiagDense(4, v)

	resQMatrix.Mul(rMat, vMat)
	resQMatrix.Mul(resQMatrix, lMat)
	sub.Sub(resQMatrix, expQMatrix)
	t.Logf("resQ = %v", mat.Formatted(resQMatrix, mat.Prefix("                 "), mat.Squeeze()))

	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			if math.Abs(sub.At(i, j)) > TEST_PRECISION {
				t.Errorf("Expected QMatrix is different from Resulting QMatrix %v", sub.At(i, j))
			}
		}
	}
}

func TestTN93Eigens(t *testing.T) {
	var v []float64
	var l, r *mat.Dense
	var vMat, lMat, rMat mat.Matrix
	var err error
	var expQMatrix *mat.Dense
	var resQMatrix *mat.Dense = mat.NewDense(4, 4, nil)
	var sub *mat.Dense = mat.NewDense(4, 4, nil)

	expQMatrix = mat.NewDense(4, 4, []float64{
		-1.0000000, 0.3333333, 0.3333333, 0.3333333,
		0.3333333, -1.0000000, 0.3333333, 0.3333333,
		0.3333333, 0.3333333, -1.0000000, 0.3333333,
		0.3333333, 0.3333333, 0.3333333, -1.0000000,
	})
	t.Logf("expQ = %v", mat.Formatted(expQMatrix, mat.Prefix("                 "), mat.Squeeze()))

	m := NewTN93Model()
	m.InitModel(1., 1., 1./4., 1./4., 1./4., 1./4.)

	if v, l, r, err = m.Eigens(); err != nil {
		t.Errorf("Error while computing F81 eigen vectors: %v", err)
	}

	// Transpose because it was in col-major formrat
	lMat = l
	rMat = r

	vMat = mat.NewDiagDense(4, v)

	resQMatrix.Mul(rMat, vMat)
	resQMatrix.Mul(resQMatrix, lMat)
	sub.Sub(resQMatrix, expQMatrix)
	t.Logf("resQ = %v", mat.Formatted(resQMatrix, mat.Prefix("                 "), mat.Squeeze()))

	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			if math.Abs(sub.At(i, j)) > TEST_PRECISION {
				t.Errorf("Expected QMatrix is different from Resulting QMatrix %v", sub.At(i, j))
			}
		}
	}
}

func TestK2PPij(t *testing.T) {

	pij := func(k float64, l float64, i, j int) float64 {
		pts := (0.25 - 0.5*math.Exp(-(2.0*k+1.0)/(k+1)*l) + 0.25*math.Exp(-2.0/(k+1.0)*l))
		ptr := 0.5 * (0.5 - 0.5*math.Exp(-2.0/(k+1)*l))
		if (i == 0 && j == 2) || (i == 1 && j == 3) ||
			(i == 2 && j == 0) || (i == 3 && j == 1) {
			// Transition
			return pts
		} else if i == j {
			return 1.0 - (pts + 2.0*ptr)
		} else {
			// Transversion
			return ptr
		}
	}

	m := NewK2PModel()
	pijmat, _ := models.NewPij(m, 1.0)

	for _, kappa := range []float64{0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0} {
		m.InitModel(kappa)
		for _, l := range []float64{0.01, 0.05, 0.1, 0.2, 0.3, 0.4} {
			pijmat.SetLength(l)
			for i := 0; i < 4; i++ {
				for j := 0; j < 4; j++ {
					if math.Abs(pijmat.Pij(i, j)-pij(0.5*kappa, l, i, j)) > TEST_PRECISION {
						t.Errorf("k=%f, l=%f | Pijmat[%d,%d] != Pijanalytial(%d,%d): %f/%f", kappa, l, i, j, i, j, pijmat.Pij(i, j), pij(kappa, l, i, j))
					}
				}
			}
		}
	}
}

func TestJCPij(t *testing.T) {

	pij := func(l float64, i, j int) float64 {
		p := 0.25 * (1.0 - math.Exp(-4.0/3.0*l))
		if i != j {
			return p
		}
		return p + math.Exp(-4.0/3.0*l)
	}

	m := NewJCModel()
	pijmat, _ := models.NewPij(m, 1.0)
	m.InitModel()
	for _, l := range []float64{0.01, 0.05, 0.1, 0.2, 0.3, 0.4} {
		pijmat.SetLength(l)
		for i := 0; i < 4; i++ {
			for j := 0; j < 4; j++ {
				if math.Abs(pijmat.Pij(i, j)-pij(l, i, j)) > TEST_PRECISION {
					t.Errorf("l=%f | Pijmat[%d,%d] != Pijanalytial(%d,%d): %f/%f", l, i, j, i, j, pijmat.Pij(i, j), pij(l, i, j))
				}
			}
		}
	}
}
