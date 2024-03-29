package protein

import (
	"fmt"
	"math"

	"gonum.org/v1/gonum/mat"
)

const (
	MODEL_DAYHOFF = iota
	MODEL_JTT
	MODEL_MTREV
	MODEL_LG
	MODEL_WAG
	MODEL_HIVB
	MODEL_AB

	BL_MIN  = 1.e-08
	BL_MAX  = 100.0
	DBL_MIN = 2.2250738585072014e-308
)

var DBL_EPSILON float64 = math.Nextafter(1, 2) - 1

type ProtModel struct {
	pi         []float64  // aa frequency
	mat        *mat.Dense // substitution matrix
	mr         float64    //MeanRate
	eigen      *mat.Eigen // Eigen Values/vectors
	leigenvect *mat.Dense // Left Eigen Vector (Inv of Eigen Vector)
	reigenvect *mat.Dense // Right Eigen Vector
	eval       []float64  // Eigen values
	alpha      float64    // Alpha
	usegamma   bool
}

// Initialize a new protein model, given the name of the model as const int:
// MODEL_DAYHOFF, MODEL_JTT, MODEL_MTREV, MODEL_LG, MODEL_WAG or MODEL_HIVB
func NewProtModel(model int, usegamma bool, alpha float64) (*ProtModel, error) {
	var m *mat.Dense
	var pi []float64
	switch model {
	case MODEL_DAYHOFF:
		m, pi = DayoffMats()
	case MODEL_JTT:
		m, pi = JTTMats()
	case MODEL_MTREV:
		m, pi = MtREVMats()
	case MODEL_LG:
		m, pi = LGMats()
	case MODEL_WAG:
		m, pi = WAGMats()
	case MODEL_HIVB:
		m, pi = HIVBMats()
	case MODEL_AB:
		m, pi = ABMats()

	default:
		return nil, fmt.Errorf("this protein model is not implemented")
	}
	return &ProtModel{
		pi,
		m,
		-1.0,
		nil,
		nil,
		nil,
		nil,
		alpha,
		usegamma,
	}, nil
}

// Returns code of the model
// If the model does not exist, returns -1
func ModelStringToInt(model string) int {
	switch model {
	case "dayoff":
		return MODEL_DAYHOFF
	case "jtt":
		return MODEL_JTT
	case "mtrev":
		return MODEL_MTREV
	case "lg":
		return MODEL_LG
	case "wag":
		return MODEL_WAG
	case "hivb":
		return MODEL_HIVB
	case "ab":
		return MODEL_AB
	default:
		return -1
	}
}

// Initialize model with given aa frequencies
// If aafreqs is nil, then uses model frequencies
func (model *ProtModel) InitModel(aafreqs []float64) error {
	var i, j int
	var sum float64
	var ok bool
	var u mat.CDense

	ns := 20

	if model.mat == nil || model.pi == nil {
		return fmt.Errorf("matrices have not been initialized")
	}

	if aafreqs != nil && len(aafreqs) != ns {
		return fmt.Errorf("aa frequency array does not have a length of 20")
	}
	if aafreqs != nil {
		model.pi = aafreqs
	}

	/* multiply the nth col of Q by the nth term of pi/100 just as in PAML */
	model.mat.Apply(func(i, j int, v float64) float64 { return v * model.pi[j] / 100.0 }, model.mat)

	/* compute diagonal terms of Q and mean rate mr = l/t */
	model.mr = .0
	for i = 0; i < ns; i++ {
		sum = .0
		for j = 0; j < ns; j++ {
			sum += model.mat.At(i, j)
		}
		model.mat.Set(i, i, -sum)
		model.mr += model.pi[i] * sum
	}

	/* scale instantaneous rate matrix so that mu=1 */
	model.mat.Apply(func(i, j int, v float64) float64 { return v / model.mr }, model.mat)

	model.eigen = &mat.Eigen{}
	if ok = model.eigen.Factorize(model.mat, mat.EigenRight); !ok {
		return fmt.Errorf("problem during matrix decomposition")
	}
	model.reigenvect = mat.NewDense(ns, ns, nil)
	model.leigenvect = mat.NewDense(20, 20, nil)
	model.eigen.VectorsTo(&u)
	model.eval = make([]float64, ns)
	for i, b := range model.eigen.Values(nil) {
		model.eval[i] = real(b)
	}
	model.reigenvect.Apply(func(i, j int, val float64) float64 { return real(u.At(i, j)) }, model.reigenvect)
	model.leigenvect.Inverse(model.reigenvect)

	return nil
}

func (model *ProtModel) Eigens() (val []float64, leftvectors, rightvectors *mat.Dense, err error) {
	return model.eval, model.leigenvect, model.reigenvect, nil
}

func (model *ProtModel) ReigenVects() (rightvectors *mat.Dense) {
	return model.reigenvect
}

func (model *ProtModel) LeigenVects() (leftvectors *mat.Dense) {
	return model.leigenvect
}

func (model *ProtModel) Eval() (val []float64) {
	return model.eval
}

func (model *ProtModel) Alpha() float64 {
	return model.alpha
}

func (model *ProtModel) UseGamma() bool {
	return model.usegamma
}

func (model *ProtModel) Pi(i int) float64 {
	return model.pi[i]
}

func (model *ProtModel) Analytical() bool {
	return false
}

func (model *ProtModel) Pij(i, j int, l float64) float64 {
	return -1.0
}

func (model *ProtModel) NState() int {
	return 20
}

func (model *ProtModel) PrintMat() {
	fa := mat.Formatted(model.mat, mat.Prefix("    "), mat.Squeeze())
	fmt.Printf("m = %v\n\n", fa)
}

func (model *ProtModel) PrintFreqs() {
	fmt.Printf("pi = [")
	for i, v := range model.pi {
		if i > 0 {
			fmt.Printf(",")
		}
		fmt.Printf("%f", v)
	}
	fmt.Printf("]\n")
}
