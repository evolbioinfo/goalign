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
	model        *protein.ProtModel
	globalAAFreq bool // Global amino acid frequency: If true we use model frequencies, else data frequencies
	removegaps   bool
	stepsize     int
}

// Initialize a new protein model, given the name of the model as const int:
// MODEL_DAYHOFF, MODEL_JTT, MODEL_MTREV, MODEL_LG or MODEL_WAG
func NewProtDistModel(model int, globalAAFreq bool, usegamma bool, alpha float64, removegaps bool) (*ProtDistModel, error) {
	m, err := protein.NewProtModel(model, usegamma, alpha)
	if err != nil {
		return nil, err
	}
	return &ProtDistModel{
		m,
		globalAAFreq,
		removegaps,
		1,
	}, nil
}

func (model *ProtDistModel) InitModel(a align.Alignment, weights []float64) (err error) {
	var pi []float64

	ns := 20

	// Count equilibrium frequencies from input alignment (do not use model frequencies)
	if !model.globalAAFreq {
		if ns = len(a.AlphabetCharacters()); ns != 20 {
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
	return model.model.Ns()
}

func (model *ProtDistModel) pMat(len float64) {
	model.model.PMat(len)
}
func (model *ProtDistModel) pij(i, j int) float64 {
	return model.model.Pij(i, j)
}

func (model *ProtDistModel) pi(i int) float64 {
	return model.model.Pi(i)
}
