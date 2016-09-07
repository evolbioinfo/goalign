package distance

import (
	"errors"
	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/fredericlemoine/goalign/stats"
)

type DistModel interface {
	InitModel(al align.Alignment)
	Distance(seq1 []rune, seq2 []rune, weigths []float64) float64
}

/* Returns the right model depending on the args */
func Model(modelType string) DistModel {
	var model DistModel
	switch modelType {
	case "jc":
		model = NewJCModel()
	case "k2p":
		model = NewK2PModel()
	case "pdist":
		model = NewPDistModel()
	case "f81":
		model = NewF81Model()
	default:
		io.ExitWithMessage(errors.New("This model is not implemented : " + modelType))
	}
	return model
}

/* Return a normalized vector of weights */
func BuildWeights(al align.Alignment) []float64 {
	outweights := make([]float64, al.Length())
	// Parameters of the binomial
	n := float64(al.Length())
	p := float64(1.0 / n)
	// Estimated parameters of the gamma
	// that fits the binomial
	alpha := (n * p / (1 - p))
	beta := 1 - p
	total := 0.0
	for i := 0; i < al.Length(); i++ {
		outweights[i] = stats.Gamma(alpha, beta)
		total += outweights[i]
	}
	// Normalizing the vector
	for i := 0; i < al.Length(); i++ {
		outweights[i] = outweights[i] * float64(al.Length()) / total
	}

	return outweights
}

/* Compute a matrix distance, with weights associated to each alignment positions */
/* If weights == nil, then all weights are considered 1 */
func DistMatrix(al align.Alignment, weights []float64, model DistModel) [][]float64 {
	outmatrix := make([][]float64, al.NbSequences())
	for i := 0; i < al.NbSequences(); i++ {
		outmatrix[i] = make([]float64, al.NbSequences())
	}
	for i := 0; i < al.NbSequences(); i++ {
		seq1, _ := al.GetSequenceChar(i)
		outmatrix[i][i] = 0
		for j := i + 1; j < al.NbSequences(); j++ {
			seq2, _ := al.GetSequenceChar(j)
			outmatrix[i][j] = model.Distance(seq1, seq2, weights)
			outmatrix[j][i] = outmatrix[i][j]
		}
	}
	return outmatrix
}

/* Returns true if it is a transition, false otherwize */
func isTransition(nt1 rune, nt2 rune) bool {
	return ((nt1 == 'A' && nt2 == 'G') || (nt1 == 'G' && nt2 == 'A') ||
		(nt1 == 'T' && nt2 == 'C') || (nt1 == 'C' && nt2 == 'T'))
}

/* Returns true if it is a transversion, false otherwize */
func isTransversion(nt1 rune, nt2 rune) bool {
	return ((nt1 == 'A' && nt2 == 'C') || (nt1 == 'C' && nt2 == 'A') ||
		(nt1 == 'G' && nt2 == 'T') || (nt1 == 'T' && nt2 == 'G') ||
		(nt1 == 'T' && nt2 == 'A') || (nt1 == 'A' && nt2 == 'T') ||
		(nt1 == 'C' && nt2 == 'G') || (nt1 == 'G' && nt2 == 'C'))
}

/* Count number of mutations and associate a weight to them */
func countMutations(seq1 []rune, seq2 []rune, weights []float64) (transitions, transversions float64) {
	transitions, transversions = 0.0, 0.0
	for i := 0; i < len(seq1); i++ {
		w := 1.0
		if weights != nil {
			w = weights[i]
		}
		if seq1[i] != seq2[i] {
			if isTransversion(seq1[i], seq2[i]) {
				transversions += float64(w)
			} else if isTransition(seq1[i], seq2[i]) {
				transitions += float64(w)
			}
		}
	}
	return
}

/* Count number of mutations and associate a weight to them */
func countDiffs(seq1 []rune, seq2 []rune, weights []float64) (nbdiffs float64) {
	nbdiffs = 0
	for i := 0; i < len(seq1); i++ {
		w := 1.0
		if weights != nil {
			w = weights[i]
		}
		if seq1[i] != seq2[i] {
			nbdiffs += float64(w)
		}
	}
	return
}
