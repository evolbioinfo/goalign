package distance

import (
	"errors"
	"fmt"
	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/fredericlemoine/goalign/stats"
	"sync"
	"unicode"
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

func init2DFloat(dim1, dim2 int) [][]float64 {
	out := make([][]float64, dim1)
	for d := 0; d < dim1; d++ {
		out[d] = make([]float64, dim2)
	}
	return out
}

/* Returns the right model depending on the args */
func Model(modelType string, removegaps bool) DistModel {
	var model DistModel
	switch modelType {
	case "jc":
		model = NewJCModel(removegaps)
	case "k2p":
		model = NewK2PModel(removegaps)
	case "pdist":
		model = NewPDistModel(removegaps)
	case "f81":
		model = NewF81Model(removegaps)
	case "tn82":
		model = NewTN82Model(removegaps)
	case "tn93":
		model = NewTN93Model(removegaps)
	case "f84":
		model = NewF84Model(removegaps)
	default:
		io.ExitWithMessage(errors.New("This model is not implemented : " + modelType))
	}
	return model
}

/* Return a normalized vector of weights following a Gamma distribution*/
func BuildWeightsGamma(al align.Alignment) []float64 {
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

/* Returns a vector of weights following a Dirichlet distribution D(n ; 1,...,1)
   with n alignment length
*/
func BuildWeightsDirichlet(al align.Alignment) []float64 {
	alpha := make([]float64, al.Length(), al.Length())
	for i := 0; i < al.Length(); i++ {
		alpha[i] = 1
	}
	outweights := stats.Dirichlet(float64(al.Length()), alpha...)
	return outweights
}

/* Compute a matrix distance, with weights associated to each alignment positions */
/* If weights == nil, then all weights are considered 1 */
func DistMatrix(al align.Alignment, weights []float64, model DistModel, cpus int) [][]float64 {
	if al.Alphabet() != align.NUCLEOTIDS {
		io.ExitWithMessage(errors.New("The alignment is not nucleotidic"))
	}
	model.InitModel(al, weights)
	distchan := make(chan seqpairdist, 100)

	outmatrix := make([][]float64, al.NbSequences())
	for i := 0; i < al.NbSequences(); i++ {
		outmatrix[i] = make([]float64, al.NbSequences())
	}

	go func() {
		for i := 0; i < al.NbSequences(); i++ {
			seq1, _ := al.GetSequenceCharById(i)
			for j := i + 1; j < al.NbSequences(); j++ {
				seq2, _ := al.GetSequenceCharById(j)
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
	n1, n2 := unicode.ToUpper(nt1), unicode.ToUpper(nt2)
	return ((n1 == 'A' && n2 == 'G') || (n1 == 'G' && n2 == 'A') ||
		(n1 == 'T' && n2 == 'C') || (n1 == 'C' && n2 == 'T'))
}

/* Returns true if it is a A<->G  */
func isAG(nt1 rune, nt2 rune) bool {
	n1, n2 := unicode.ToUpper(nt1), unicode.ToUpper(nt2)
	return ((n1 == 'A' && n2 == 'G') || (n1 == 'G' && n2 == 'A'))
}

/* Returns true if it is a A<->G  */
func isCT(nt1 rune, nt2 rune) bool {
	n1, n2 := unicode.ToUpper(nt1), unicode.ToUpper(nt2)
	return ((n1 == 'T' && n2 == 'C') || (n1 == 'C' && n2 == 'T'))
}

/* Returns true if it is a transversion, false otherwize */
func isTransversion(nt1 rune, nt2 rune) bool {
	n1, n2 := unicode.ToUpper(nt1), unicode.ToUpper(nt2)
	return ((n1 == 'A' && n2 == 'C') || (n1 == 'C' && n2 == 'A') ||
		(n1 == 'G' && n2 == 'T') || (n1 == 'T' && n2 == 'G') ||
		(n1 == 'T' && n2 == 'A') || (n1 == 'A' && n2 == 'T') ||
		(n1 == 'C' && n2 == 'G') || (n1 == 'G' && n2 == 'C'))
}

/* Count number of mutations and associate a weight to them */
func countMutations(seq1 []rune, seq2 []rune, selectedSites []bool, weights []float64) (transitions, transversions, ag, ct float64, total float64) {
	transitions, transversions = 0.0, 0.0
	total = 0.0
	ag = 0.0
	ct = 0.0
	for i := 0; i < len(seq1); i++ {
		w := 1.0
		if weights != nil {
			w = weights[i]
		}
		if isNuc(seq1[i]) && isNuc(seq2[i]) && selectedSites[i] {
			if seq1[i] != seq2[i] {
				if isTransversion(seq1[i], seq2[i]) {
					transversions += w
				} else if isTransition(seq1[i], seq2[i]) {
					transitions += w
				}
				if isAG(seq1[i], seq2[i]) {
					ag += w
				} else if isCT(seq1[i], seq2[i]) {
					ct += w
				}
			}
			total += w
		}
	}
	return
}

/* Count number of mutations and associate a weight to them */
func countDiffs(seq1 []rune, seq2 []rune, selectedSites []bool, weights []float64) (nbdiffs float64, total float64) {
	nbdiffs = 0.0
	total = 0.0
	for i := 0; i < len(seq1); i++ {
		w := 1.0
		if weights != nil {
			w = weights[i]
		}
		if isNuc(seq1[i]) && isNuc(seq2[i]) && selectedSites[i] {
			if seq1[i] != seq2[i] {
				nbdiffs += float64(w)
			}
			total += w
		}
	}
	return
}

/*
Returns the proba of each nts
A=0
C=1
G=2
T=3
*/
func probaNt(al align.Alignment, selectedSites []bool, weights []float64) []float64 {
	pi := make([]float64, 4)
	total := 0.0

	w := 1.0
	for pos := 0; pos < al.Length(); pos++ {
		if weights != nil {
			w = weights[pos]
		}
		for i := 0; i < al.NbSequences(); i++ {
			seq1, _ := al.GetSequenceCharById(i)
			char := seq1[pos]
			if selectedSites[pos] {
				if isNuc(seq1[pos]) {
					pi[indexNt(char)] += w
				}
				total += w
			}
		}
	}

	for i, _ := range pi {
		pi[i] /= total
	}
	return pi
}

/*Returns the proba of each nts for the 2 sequences considered
A=0
C=1
G=2
T=3
*/
func probaNt2Seqs(seq1 []rune, seq2 []rune, selectedSites []bool, weights []float64) []float64 {
	pi := make([]float64, 4)
	total := 0.0

	w := 1.0
	for pos := 0; pos < len(seq1); pos++ {
		if weights != nil {
			w = weights[pos]
		}
		if selectedSites[pos] {
			if isNuc(seq1[pos]) && isNuc(seq2[pos]) {
				pi[indexNt(seq1[pos])] += w
				pi[indexNt(seq2[pos])] += w
			}
			total += 2 * w
		}
	}

	for i, _ := range pi {
		pi[i] /= total
	}
	return pi
}

/* Compute freq (weighted) of all pairs of nt in all pairs of sequences */
func probaNtPairs(al align.Alignment, selectedSites []bool, weights []float64) [][]float64 {
	psi := init2DFloat(4, 4)
	total := 0.0
	for i := 0; i < al.NbSequences(); i++ {
		seq1, _ := al.GetSequenceCharById(i)
		for j := i + 1; j < al.NbSequences(); j++ {
			seq2, _ := al.GetSequenceCharById(j)
			total += countNtPairs2Seq(seq1, seq2, selectedSites, weights, psi)
		}
	}

	for i, _ := range psi {
		for j, _ := range psi[i] {
			psi[i][j] /= total
		}
	}
	return psi
}

/* Compute freq (weighted) of all pairs of nt in this pair of sequences */
func countNtPairs2Seq(seq1, seq2 []rune, selectedSites []bool, weights []float64, psi [][]float64) float64 {
	total := 0.0
	w := 1.0
	for pos, char1 := range seq1 {
		if weights != nil {
			w = weights[pos]
		}
		if selectedSites[pos] {
			if isNuc(char1) && isNuc(seq2[pos]) {
				id1 := indexNt(char1)
				id2 := indexNt(seq2[pos])
				psi[id1][id2] += w
				psi[id2][id1] += w
			}
			total += w
		}
	}
	return total
}

/*
Returns the index of each nts
A=0
C=1
G=2
T=3
*/
func indexNt(nt rune) int {
	switch unicode.ToUpper(nt) {
	case 'A':
		return 0
	case 'C':
		return 1
	case 'G':
		return 2
	case 'T':
		return 3
	}
	io.ExitWithMessage(errors.New(fmt.Sprintf("No index for character: %c", nt)))
	return -1
}

func isNuc(r rune) bool {
	up := unicode.ToUpper(r)
	return up == 'A' || up == 'C' || up == 'G' || up == 'T'
}

/* Returns the sites of the alignments that contains only nucleotides and no gaps */
func selectedSites(al align.Alignment, weights []float64, removeGappedPositions bool) (float64, []bool) {
	selectedSites := make([]bool, al.Length())
	numSites := 0.0
	for l := 0; l < al.Length(); l++ {
		w := 1.0
		if weights != nil {
			w = weights[l]
		}
		selectedSites[l] = true
		for i := 0; i < al.NbSequences() && removeGappedPositions; i++ {
			seq, _ := al.GetSequenceCharById(i)
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
