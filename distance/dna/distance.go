package dna

import (
	"errors"
	"math"
	"sync"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/stats"
)

const (
	NT_DIST_OVER = 100000
)

type DistModel interface {
	InitModel(al align.Alignment, weights []float64, gamma bool, alpha float64) error
	Distance(seq1 []int, seq2 []int, weigths []float64) (float64, error)
	Sequence(i int) ([]int, error)
}

type seqpairdist struct {
	i, j       int
	seq1, seq2 []int // Sequences encoded with align.Nt2IndexIUPAC
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
func Model(modelType string, removegaps bool) (model DistModel, err error) {
	switch modelType {
	case "jc":
		model = NewJCModel(removegaps)
	case "k2p":
		model = NewK2PModel(removegaps)
	case "pdist":
		model = NewPDistModel(removegaps)
	case "rawdist":
		model = NewRawDistModel(removegaps)
	case "f81":
		model = NewF81Model(removegaps)
	// case "tn82":
	// 	model = NewTN82Model(removegaps)
	case "tn93":
		model = NewTN93Model(removegaps)
	case "f84":
		model = NewF84Model(removegaps)
	default:
		err = errors.New("This model is not implemented : " + modelType)
	}
	return
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
func DistMatrix(al align.Alignment, weights []float64, model DistModel, gamma bool, alpha float64, cpus int) (outmatrix [][]float64, err error) {
	if al.Alphabet() != align.NUCLEOTIDS {
		err = errors.New("The alignment is not nucleotidic")
		return
	}
	if err = model.InitModel(al, weights, gamma, alpha); err != nil {
		return
	}
	distchan := make(chan seqpairdist, 100)

	// Chan of comparisons that gives infty or <0 distance
	uncompute := make([]seqpairdist, 0, 100)
	var mux sync.Mutex

	outmatrix = make([][]float64, al.NbSequences())
	for i := 0; i < al.NbSequences(); i++ {
		outmatrix[i] = make([]float64, al.NbSequences())
	}

	go func() {
		defer close(distchan)
		var seq1, seq2 []int
		for i := 0; i < al.NbSequences(); i++ {
			if seq1, err = model.Sequence(i); err != nil {
				return
			}
			for j := i + 1; j < al.NbSequences(); j++ {
				if seq2, err = model.Sequence(j); err != nil {
					return
				}
				distchan <- seqpairdist{i, j, seq1, seq2, model, weights}
			}
		}
	}()
	if err != nil {
		return
	}

	var wg sync.WaitGroup
	max := 0.0
	for cpu := 0; cpu < cpus; cpu++ {
		wg.Add(1)
		go func() {
			for sp := range distchan {
				if sp.i == sp.j {
					outmatrix[sp.i][sp.i] = 0
				} else {
					if outmatrix[sp.i][sp.j], err = model.Distance(sp.seq1, sp.seq2, sp.weights); err != nil {
						return
					}
					outmatrix[sp.j][sp.i] = outmatrix[sp.i][sp.j]
					mux.Lock()
					if outmatrix[sp.i][sp.j] < 0 || outmatrix[sp.i][sp.j] == math.Inf(1) || outmatrix[sp.i][sp.j] > NT_DIST_OVER {
						uncompute = append(uncompute, seqpairdist{sp.i, sp.j, nil, nil, nil, nil})
					} else if outmatrix[sp.i][sp.j] > max {
						max = outmatrix[sp.i][sp.j]
					}
					mux.Unlock()
				}
			}
			wg.Done()
		}()
	}
	wg.Wait()

	for _, sp := range uncompute {
		outmatrix[sp.i][sp.j] = 2 * max
		outmatrix[sp.j][sp.i] = 2 * max
	}

	return
}

/* Returns true if it is a transition, false otherwize */
func isTransition(n1 int, n2 int) bool {
	return ((n1 == align.NT_A && n2 == align.NT_G) || (n1 == align.NT_G && n2 == align.NT_A) ||
		(n1 == align.NT_T && n2 == align.NT_C) || (n1 == align.NT_C && n2 == align.NT_T))
}

/* Returns true if it is a A<->G  */
func isAG(n1 int, n2 int) bool {
	return ((n1 == align.NT_A && n2 == align.NT_G) || (n1 == align.NT_G && n2 == align.NT_A))
}

/* Returns true if it is a C<->T  */
func isCT(n1 int, n2 int) bool {
	return ((n1 == align.NT_T && n2 == align.NT_C) || (n1 == align.NT_C && n2 == align.NT_T))
}

/* Returns true if it is a transversion, false otherwize */
func isTransversion(n1 int, n2 int) bool {
	return ((n1 == align.NT_A && n2 == align.NT_C) || (n1 == align.NT_C && n2 == align.NT_A) ||
		(n1 == align.NT_G && n2 == align.NT_T) || (n1 == align.NT_T && n2 == align.NT_G) ||
		(n1 == align.NT_T && n2 == align.NT_A) || (n1 == align.NT_A && n2 == align.NT_T) ||
		(n1 == align.NT_C && n2 == align.NT_G) || (n1 == align.NT_G && n2 == align.NT_C) ||
		(n1 == align.NT_R && n2 == align.NT_Y) || (n1 == align.NT_Y && n2 == align.NT_R))
}

/* Count number of mutations and associate a weight to them */
func countMutations(seq1 []int, seq2 []int, selectedSites []bool, weights []float64) (transitions, transversions, ag, ct float64, total float64) {
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
func countDiffs(seq1 []int, seq2 []int, selectedSites []bool, weights []float64) (nbdiffs float64, total float64) {
	nbdiffs = 0.0
	total = 0.0
	for i := 0; i < len(seq1); i++ {
		w := 1.0
		if weights != nil {
			w = weights[i]
		}
		if isNuc(seq1[i]) && isNuc(seq2[i]) && selectedSites[i] {
			if seq1[i] != seq2[i] {
				diff, _ := align.NtIUPACDifference(seq1[i], seq2[i])
				diff = diff * w
				nbdiffs += diff
			}
			total += w
		}
	}
	return
}

/* Count number of mutations (including gaps to nt) */
func countDiffsWithGaps(seq1, seq2 []int, selectedSites []bool, weights []float64) (nbdiffs float64, total float64) {
	nbdiffs = 0.0
	total = 0.0
	for i := 0; i < len(seq1); i++ {
		w := 1.0
		if weights != nil {
			w = weights[i]
		}

		if (isNuc(seq1[i]) || isNuc(seq2[i])) && selectedSites[i] {
			if seq1[i] != seq2[i] {
				diff, _ := align.NtIUPACDifference(seq1[i], seq2[i])
				diff = diff * w
				nbdiffs += diff
			}
			total += w
		}
	}
	return
}

/* Count number of mutations (including gaps to nt when they are internal, not on the border) */
func countDiffsWithInternalGaps(seq1, seq2 []int, selectedSites []bool, weights []float64) (nbdiffs float64, total float64) {
	nbdiffs = 0.0
	total = 0.0
	firstgaps1 := true
	firstgaps2 := true
	tmpgapdiffs1 := 0.0
	tmpgapdiffs2 := 0.0

	for i := 0; i < len(seq1); i++ {
		w := 1.0
		if weights != nil {
			w = weights[i]
		}
		firstgaps1 = firstgaps1 && !isNuc(seq1[i])
		firstgaps2 = firstgaps2 && !isNuc(seq2[i])

		if (isNuc(seq1[i]) || isNuc(seq2[i])) && (!firstgaps1 && !firstgaps2) {
			if seq1[i] != seq2[i] {
				diff, _ := align.NtIUPACDifference(seq1[i], seq2[i])
				diff = diff * w
				nbdiffs += diff
				tmpgapdiffs1 += diff
				tmpgapdiffs2 += diff
			}
			total += w

			if isNuc(seq1[i]) {
				tmpgapdiffs1 = .0
			}

			if isNuc(seq2[i]) {
				tmpgapdiffs2 = .0
			}
		}
	}

	nbdiffs -= math.Max(tmpgapdiffs1, tmpgapdiffs2)
	total -= math.Max(tmpgapdiffs1, tmpgapdiffs2)

	return
}

/*
Returns the proba of each nts
A=0
C=1
G=2
T=3
*/
func probaNt(sequenceCodes [][]int, selectedSites []bool, weights []float64) ([]float64, error) {
	var idx []int
	var err error
	var i, l, seqidx, pos int
	var seq1 []int

	pi := make([]float64, 4)
	total := 0.0

	if len(sequenceCodes) > 0 {
		l = len(sequenceCodes[0])
	}

	w := 1.0
	for pos = 0; pos < l; pos++ {
		if weights != nil {
			w = weights[pos]
		}
		for seqidx = 0; seqidx < len(sequenceCodes); seqidx++ {
			seq1 = sequenceCodes[seqidx]
			if selectedSites[pos] {
				if isNuc(seq1[pos]) {
					if idx, err = align.PossibleNtIUPAC(seq1[pos]); err != nil {
						return nil, err
					}
					for _, i = range idx {
						pi[i] += w / float64(len(idx))
					}
				}
				total += w
			}
		}
	}

	for i = range pi {
		pi[i] /= total
	}
	return pi, nil
}

/*Returns the proba of each nts for the 2 sequences considered
A=0
C=1
G=2
T=3
*/
func probaNt2Seqs(seq1 []int, seq2 []int, selectedSites []bool, weights []float64) ([]float64, error) {
	var idx []int
	var i int
	var err error

	pi := make([]float64, 4)
	total := 0.0

	w := 1.0
	for pos := 0; pos < len(seq1); pos++ {
		if weights != nil {
			w = weights[pos]
		}
		if selectedSites[pos] {
			if isNuc(seq1[pos]) && isNuc(seq2[pos]) {

				if idx, err = align.PossibleNtIUPAC(seq1[pos]); err != nil {
					return nil, err
				}
				for _, i = range idx {
					pi[i] += w / float64(len(idx))
				}

				if idx, err = align.PossibleNtIUPAC(seq2[pos]); err != nil {
					return nil, err
				}
				for _, i = range idx {
					pi[i] += w / float64(len(idx))
				}
			}
			total += 2 * w
		}
	}

	for i = range pi {
		pi[i] /= total
	}
	return pi, nil
}

/* Compute freq (weighted) of all pairs of nt in all pairs of sequences */
func probaNtPairs(sequenceCodes [][]int, selectedSites []bool, weights []float64) ([][]float64, error) {
	var seq1, seq2 []int

	psi := init2DFloat(4, 4)
	total := 0.0

	for i := 0; i < len(sequenceCodes); i++ {
		seq1 = sequenceCodes[i]
		for j := i + 1; j < len(seq1); j++ {
			seq2 = sequenceCodes[j]
			if nb, err := countNtPairs2Seq(seq1, seq2, selectedSites, weights, psi); err == nil {
				total += nb
			} else {
				return nil, err
			}
		}
	}

	for i := range psi {
		for j := range psi[i] {
			psi[i][j] /= total
		}
	}
	return psi, nil
}

/* Compute freq (weighted) of all pairs of nt in this pair of sequences */
func countNtPairs2Seq(seq1, seq2 []int, selectedSites []bool, weights []float64, psi [][]float64) (total float64, err error) {
	total = 0.0
	w := 1.0
	for pos, char1 := range seq1 {
		if weights != nil {
			w = weights[pos]
		}
		if selectedSites[pos] {
			if isNucStrict(char1) && isNucStrict(seq2[pos]) {
				var id1, id2 []int

				if id1, err = align.PossibleNtIUPAC(char1); err != nil {
					return
				}
				if id2, err = align.PossibleNtIUPAC(seq2[pos]); err != nil {
					return
				}
				nb := float64(len(id1) * len(id2))
				for _, i1 := range id1 {
					for _, i2 := range id2 {
						psi[i1][i2] += w / nb
						if i1 != i2 {
							psi[i2][i1] += w / nb
						}
					}
				}
			}
			total += w
		}
	}
	return
}

func isNuc(r int) bool {
	return r >= 0 && r < len(align.IupacCode)
}

func isNucStrict(r int) bool {
	return r >= 0 && r < 4
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
			if al.AlphabetCharToIndex(seq[l]) == -1 || seq[l] == '*' || seq[l] == '?' || seq[l] == '-' {
				selectedSites[l] = false
				break
			}
		}
		if selectedSites[l] {
			numSites += w
		}
	}
	return numSites, selectedSites
}

func alignmentToCodes(al align.Alignment) (sequencesInCode [][]int, err error) {
	var i int
	// Sequences coded in NT_A-NT_OTHER
	sequencesInCode = make([][]int, al.NbSequences())
	i = 0
	al.IterateChar(func(name string, seq []rune) bool {
		sequencesInCode[i] = make([]int, al.Length())
		for l, r := range seq {
			if sequencesInCode[i][l], err = align.Nt2IndexIUPAC(r); err != nil {
				return true
			}
		}
		i++
		return false
	})
	return
}
