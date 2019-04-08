package protein

import (
	"fmt"
	"math"

	"github.com/gonum/matrix/mat64"
)

const (
	MODEL_DAYHOFF = iota
	MODEL_JTT
	MODEL_MTREV
	MODEL_LG
	MODEL_WAG

	PROT_DIST_MAX = 20.00
	BL_MIN        = 1.e-08
	BL_MAX        = 100.0
)

var DBL_EPSILON float64 = math.Nextafter(1, 2) - 1

type ProtModel struct {
	pi            []float64    // aa frequency
	mat           *mat64.Dense // substitution matrix
	globalAAFreq  bool         // Global amino acid frequency: If true we use model frequencies, else data frequencies
	mr            float64      //MeanRate
	eigen         *mat64.Eigen // Eigen Values/vectors
	leigenvect    *mat64.Dense // Left Eigen Vector (Inv of Eigen Vector)
	ns            int          // Number of states in the model
	pij           *mat64.Dense // Matrix of Pij
	alpha         float64      // Alpha
	stepsize      int
	n_categ       int // gamma categories // TJS 1...
	gamma_rr      float64
	gamma_r_proba float64
	usegamma      bool
}

// Initialize a new protein model, given the name of the model as const int:
// MODEL_DAYHOFF, MODEL_JTT, MODEL_MTREV, MODEL_LG or MODEL_WAG
func NewProtModel(model int, globalAAFreq bool, usegamma bool, alpha float64) (*ProtModel, error) {
	var mat *mat64.Dense
	var pi []float64
	switch model {
	case MODEL_DAYHOFF:
		mat, pi = DayoffMats()
	case MODEL_JTT:
		mat, pi = JTTMats()
	case MODEL_MTREV:
		mat, pi = MtREVMats()
	case MODEL_LG:
		mat, pi = LGMats()
	case MODEL_WAG:
		mat, pi = WAGMats()
	default:
		return nil, fmt.Errorf("This protein model is not implemented")
	}
	return &ProtModel{
		pi,
		mat,
		globalAAFreq,
		-1.0,
		nil,
		nil,
		len(pi),
		mat64.NewDense(len(pi), len(pi), nil),
		alpha,
		1,
		1,
		.0,
		.0,
		usegamma,
	}, nil
}