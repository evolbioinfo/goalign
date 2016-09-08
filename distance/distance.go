package distance

import (
	"errors"
	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/fredericlemoine/goalign/stats"
	"sync"
)

type DistModel interface {
	InitModel(al align.Alignment, weights []float64)
	Distance(seq1 []rune, seq2 []rune, weigths []float64) float64
}

type seqpairdist struct {
	i, j       int
	seq1, seq2 []rune
	model      DistModel
	weights    []float64
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
func DistMatrix(al align.Alignment, weights []float64, model DistModel, cpus int) [][]float64 {
	model.InitModel(al, weights)
	distchan := make(chan seqpairdist, 100)

	outmatrix := make([][]float64, al.NbSequences())
	for i := 0; i < al.NbSequences(); i++ {
		outmatrix[i] = make([]float64, al.NbSequences())
	}

	go func() {
		for i := 0; i < al.NbSequences(); i++ {
			seq1, _ := al.GetSequenceChar(i)
			for j := i + 1; j < al.NbSequences(); j++ {
				seq2, _ := al.GetSequenceChar(j)
				distchan <- seqpairdist{i, j, seq1, seq2, model, weights}
			}
		}
		close(distchan)
	}()

	var wg sync.WaitGroup
	for cpu := 0; cpu < cpus; cpu++ {
		wg.Add(1)
		go func() {
			for sp := range distchan {
				if sp.i == sp.j {
					outmatrix[sp.i][sp.i] = 0
				} else {
					outmatrix[sp.i][sp.j] = model.Distance(sp.seq1, sp.seq2, sp.weights)
					outmatrix[sp.j][sp.i] = outmatrix[sp.i][sp.j]
				}
			}
			wg.Done()
		}()
	}
	wg.Wait()

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
func countMutations(seq1 []rune, seq2 []rune, weights []float64) (transitions, transversions float64, total float64) {
	transitions, transversions = 0.0, 0.0
	total = 0.0
	for i := 0; i < len(seq1); i++ {
		w := 1.0
		if weights != nil {
			w = weights[i]
		}
		if isNuc(seq1[i]) && isNuc(seq2[i]) {
			if seq1[i] != seq2[i] {
				if isTransversion(seq1[i], seq2[i]) {
					transversions += float64(w)
				} else if isTransition(seq1[i], seq2[i]) {
					transitions += float64(w)
				}
			}
			total += w
		}
	}
	return
}

/* Count number of mutations and associate a weight to them */
func countDiffs(seq1 []rune, seq2 []rune, weights []float64) (nbdiffs float64, total float64) {
	nbdiffs = 0
	total = 0
	for i := 0; i < len(seq1); i++ {
		w := 1.0
		if weights != nil {
			w = weights[i]
		}
		if isNuc(seq1[i]) && isNuc(seq2[i]) {
			if seq1[i] != seq2[i] {
				nbdiffs += float64(w)
			}
			total += w
		}
	}
	return
}

func isNuc(r rune) bool {
	return r == 'A' || r == 'C' || r == 'G' || r == 'T'
}

/* Returns the sites of the alignments that contains only nucleotides and no gaps */
func selectedSites(al align.Alignment, weights []float64) (float64, []bool) {
	selectedSites := make([]bool, al.Length())
	numSites := 0.0
	for l := 0; l < al.Length(); l++ {
		w := 1.0
		if weights != nil {
			w = weights[l]
		}
		selectedSites[l] = true
		for i := 0; i < al.NbSequences(); i++ {
			seq, _ := al.GetSequenceChar(i)
			if !isNuc(seq[l]) || seq[l] == '*' || seq[l] == '?' || seq[l] == '-' {
				selectedSites[l] = false
			}
		}
		if selectedSites[l] {
			numSites += w
		}
	}
	return numSites, selectedSites
}
