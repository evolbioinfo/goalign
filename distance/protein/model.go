package protein

import (
	"fmt"
	"math"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/gotree/models/protein"
	"gonum.org/v1/gonum/mat"
)

const (
	PROT_DIST_MAX = 20.00
)

type ProtDistModel struct {
	model      *protein.ProtModel
	modelfreqs bool //  amino acid frequency: If true we use model frequencies, else empirical data frequencies
	removegaps bool
	pij        *mat.Dense // Matrix of Pij
	stepsize   int
}

// Initialize a new protein model, given the name of the model as const int:
// MODEL_DAYHOFF, MODEL_JTT, MODEL_MTREV, MODEL_LG or MODEL_WAG
func NewProtDistModel(model int, modelfreqs bool, usegamma bool, alpha float64, removegaps bool) (*ProtDistModel, error) {
	m, err := protein.NewProtModel(model, usegamma, alpha)
	if err != nil {
		return nil, err
	}
	return &ProtDistModel{
		m,
		modelfreqs,
		removegaps,
		nil,
		1,
	}, nil
}

func (model *ProtDistModel) InitModel(a align.Alignment, weights []float64) (err error) {
	var pi []float64

	ns := model.Ns()

	// Count equilibrium frequencies from input alignment (do not use model frequencies)
	if !model.modelfreqs {
		if ns = len(a.AlphabetCharacters()); ns != model.Ns() {
			err = fmt.Errorf("Alphabet has not a length of 20")
			return
		}
		if weights == nil {
			weights = make([]float64, a.Length())
			for i, _ := range weights {
				weights[i] = 1.0
			}
		}
		_, selected := selectedSites(a, weights, model.removegaps)
		if pi, err = aaFrequency(a, weights, selected); err != nil {
			return
		}
	}
	model.pij = mat.NewDense(len(pi), len(pi), nil)
	model.model.InitModel(pi)
	return nil
}

// Basic JC69 Protein Distance Matrix
func (model *ProtDistModel) JC69Dist(a align.Alignment, weights []float64, selected []bool) (p *mat.Dense, q *mat.Dense, dist *mat.Dense) {
	var site, i, j, k int
	var len *mat.Dense
	ns := model.Ns()

	len = mat.NewDense(a.NbSequences(), a.NbSequences(), nil)
	p = mat.NewDense(a.NbSequences(), a.NbSequences(), nil)
	q = mat.NewDense(a.NbSequences(), a.NbSequences(), nil)
	dist = mat.NewDense(a.NbSequences(), a.NbSequences(), nil)

	for site = 0; site < a.Length(); site += model.stepsize {
		if selected[site] {
			for j = 0; j < a.NbSequences()-1; j++ {
				s1, _ := a.GetSequenceCharById(j)
				for k = j + 1; k < a.NbSequences(); k++ {
					s2, _ := a.GetSequenceCharById(k)
					if (!isAmbigu(s1[site])) && (!isAmbigu(s2[site])) {
						len.Set(j, k, len.At(j, k)+weights[site])
						len.Set(k, j, weights[site])
						for n, c1 := range s1[site : site+model.stepsize] {
							if c1 != s2[site+n] {
								p.Set(j, k, p.At(j, k)+weights[site])
								break
							}
						}
					}
				}
			}
		}
	}

	for i = 0; i < a.NbSequences()-1; i++ {
		for j = i + 1; j < a.NbSequences(); j++ {
			if len.At(i, j) > 0 {
				p.Set(i, j, p.At(i, j)/len.At(i, j))
			} else {
				p.Set(i, j, 1.)
			}

			p.Set(j, i, p.At(i, j))

			if (1. - float64(ns)/float64(ns-1.)*p.At(i, j)) < .0 {
				dist.Set(i, j, PROT_DIST_MAX)
			} else {
				dist.Set(i, j, -float64(ns-1.)/float64(ns)*math.Log(1.-float64(ns)/float64(ns-1.)*p.At(i, j)))
			}
			if dist.At(i, j) > PROT_DIST_MAX {
				dist.Set(i, j, PROT_DIST_MAX)
			}
			dist.Set(j, i, dist.At(i, j))
		}
	}

	len = nil

	return p, q, dist
}

func (model *ProtDistModel) Ns() int {
	return 20
}

func (model *ProtDistModel) pMat(l float64) {
	if l < BL_MIN {
		model.pMatZeroBrLen()
	} else {
		model.pMatEmpirical(l)
	}
}

func (model *ProtDistModel) pMatZeroBrLen() {
	model.pij.Apply(func(i, j int, v float64) float64 {
		if i == j {
			return 1.0
		}
		return 0.0
	}, model.pij)
}

/********************************************************************/
/*                    Code taken from FastME                        */
/* Computes the substitution probability matrix
 * from the initial substitution rate matrix and frequency vector
 * and one specific branch length
 *
 * input : l , branch length
 * input : mod , choosen model parameters, qmat and pi
 * ouput : Pij , substitution probability matrix
 *
 * matrix P(l) is computed as follows :
 * P(l) = exp(Q*t) , where :
 *
 *   Q = substitution rate matrix = Vr*D*inverse(Vr) , where :
 *
 *     Vr = right eigenvector matrix for Q
 *     D  = diagonal matrix of eigenvalues for Q
 *
 *   t = time interval = l / mr , where :
 *
 *     mr = mean rate = branch length/time interval
 *        = sum(i)(pi[i]*p(i->j)) , where :
 *
 *       pi = state frequency vector
 *       p(i->j) = subst. probability from i to a different state
 *               = -Q[ii] , as sum(j)(Q[ij]) +Q[ii] = 0
 *
 * the Taylor development of exp(Q*t) gives :
 * P(l) = Vr*exp(D*t)        *inverse(Vr)
 *      = Vr*pow(exp(D/mr),l)*inverse(Vr)
 *
 * for performance we compute only once the following matrices :
 * Vr, inverse(Vr), exp(D/mr)
 * thus each time we compute P(l) we only have to :
 * make 20 times the operation pow()
 * make 2 20x20 matrix multiplications, that is :
 *   16000 = 2x20x20x20 times the operation *
 *   16000 = 2x20x20x20 times the operation +
 *   which can be reduced to (the central matrix being diagonal) :
 *   8400 = 20x20 + 20x20x20 times the operation *
 *   8000 = 20x20x20 times the operation + */
func (model *ProtDistModel) pMatEmpirical(len float64) {
	var i, k int
	var U, V *mat.Dense
	var R []float64
	var expt []float64
	var uexpt *mat.Dense
	var tmp float64

	ns := model.Ns()
	U = model.model.ReigenVects() //mod->eigen->r_e_vect;
	R = model.model.Eval()        //mod->eigen->e_val;// To take only real part from that vector /* eigen value matrix */
	V = model.model.LeigenVects()
	expt = make([]float64, ns)        //model.eigen.Values(nil) // To take only imaginary part from that vector
	uexpt = mat.NewDense(ns, ns, nil) //model.eigen.Vectors() //  don't know yet how to handle that // mod->eigen->r_e_vect_im;

	model.pij.Apply(func(i, j int, v float64) float64 { return .0 }, model.pij)
	tmp = .0

	for k = 0; k < ns; k++ {
		expt[k] = R[k]
	}

	alpha := model.model.Alpha()
	if model.model.UseGamma() && (math.Abs(alpha) > DBL_EPSILON) {
		// compute pow (alpha / (alpha - e_val[i] * l), alpha)
		for i = 0; i < ns; i++ {
			tmp = alpha / (alpha - (R[i] * len))
			expt[i] = math.Pow(tmp, alpha)
		}
	} else {
		for i = 0; i < ns; i++ {
			expt[i] = float64(math.Exp(R[i] * len))
		}
	}

	// multiply Vr* pow (alpha / (alpha - e_val[i] * l), alpha) *Vi into Pij
	uexpt.Apply(func(i, j int, v float64) float64 {
		return U.At(i, j) * expt[j]
	}, uexpt)
	model.pij.Apply(func(i, j int, v float64) float64 {
		for k = 0; k < ns; k++ {
			v += uexpt.At(i, k) * V.At(k, j)
		}
		if v < DBL_MIN {
			v = DBL_MIN
		}
		return v

	}, model.pij)
}
