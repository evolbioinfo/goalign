package align

import (
	"fmt"
)

type PairwiseAligner interface {
	AlignEnds() (int, int)
	AlignStarts() (int, int)
	SetGapOpenScore(open float64)
	SetGapExtendScore(extend float64)
	SetScore(match, mismatch float64)
	MaxScore() float64
	Alignment() (Alignment, error)
	AlignmentStr() string
}

const (
	ALIGN_UP = iota
	ALIGN_LEFT
	ALIGN_DIAG
	ALIGN_STOP

	ALIGN_ALGO_SW = iota
	ALIGN_ALGO_ATG
)

// Taken from EMBOSS water
//
// This matrix was created by Todd Lowe   12/10/92
//
// Uses ambiguous nucleotide codes, probabilities rounded to
//  nearest integer
//
// Lowest score = -4, Highest score = 5
//

type pwaligner struct {
	algo   int
	matrix [][]float64
	trace  [][]int   // trace matrix
	maxa   []float64 // keep track of best gap opened

	maxscore         float64 // Maximum score of the matrix
	maxi, maxj       int     // Indices of the maximum score
	start1, start2   int
	end1, end2       int
	seq1, seq2       Sequence
	seq1ali, seq2ali []rune // aligned sequences
	alistr           []rune // comparison between sequences
	gapopen          float64
	gapextend        float64
	match            float64
	mismatch         float64
	submatrix        [][]float64  // substitution matrix
	chartopos        map[rune]int // char to position in subst matrix

}

func NewPwAligner(seq1, seq2 Sequence, algo int) *pwaligner {
	var mat [][]float64
	var chartopos map[rune]int
	a1 := seq1.DetectAlphabet()
	a2 := seq2.DetectAlphabet()

	if a1 == NUCLEOTIDS && a2 == NUCLEOTIDS {
		mat = dnafull_subst_matrix
		chartopos = dna_to_matrix_pos
	}
	if a1 == AMINOACIDS && a2 == AMINOACIDS {
		mat = blosum62_subst_matrix
		chartopos = prot_to_matrix_pos
	}

	return &pwaligner{
		algo:      algo,
		matrix:    nil,
		trace:     nil,
		maxa:      nil,
		maxscore:  0.0,
		maxi:      0,
		maxj:      0,
		start1:    0,
		start2:    0,
		end1:      0,
		end2:      0,
		seq1:      seq1,
		seq2:      seq2,
		gapopen:   -10.0,
		gapextend: -0.5,
		match:     1.0,
		mismatch:  -1.0,
		submatrix: mat,
		chartopos: chartopos,
	}
}

func (a *pwaligner) initMatrix(l1, l2 int) {
	var i, j int

	a.matrix = make([][]float64, l1)
	a.trace = make([][]int, l1)
	a.maxa = make([]float64, l2)
	for i, _ = range a.matrix {
		a.matrix[i] = make([]float64, l2)
		a.trace[i] = make([]int, l2)
	}
}

func (a *pwaligner) SetGapScore(gap float64) {
	a.gap = gap
}

// Sets manually match and mismatch scores
// Substitution matrix is not used any more
func (a *pwaligner) SetScore(match, mismatch float64) {
	a.match = match
	a.mismatch = mismatch
	a.submatrix = nil
	a.chartopos = nil
}

func (a *pwaligner) MaxScore() float64 {
	return a.maxscore
}

func (a *pwaligner) fillMatrix() (err error) {
	switch a.algo {
	case ALIGN_ALGO_ATG:
		err = a.fillMatrix_ATG()
	default:
		err = a.fillMatrix_SW()
	}
	return
}

func (a *pwaligner) fillMatrix_SW() (err error) {
	var i, j int
	var c1, c2 rune
	var maxscore float64
	var trace int
	var leftscore, upscore, diagscore float64
	var indexseq1, indexseq2 []int // convert characters to subst matrix positions
	a.initMatrix(a.seq1.Length(), a.seq2.Length())

	// We convert characters to indices in subst matrices
	// once for all
	if indexseq1, err = a.seqToindices(a.seq1); err != nil {
		return
	}
	if indexseq2, err = a.seqToindices(a.seq2); err != nil {
		return
	}

	// We initialize first row and first column of the matrix
	var match, fnew float64
	// First row
	for i = 0; i <= l1; i++ {
		c1 = a.seq1.SequenceChar(i)
		c2 = a.seq2.SequenceChar(0)
		match = a.matchScore(c1, c2, indexseq1[i], indexseq2[0])
		fnew = 0.0
		if i > 0 {
			fnew = a.matrix[i-1][0]
			if a.trace[i-1][0] == ALIGN_LEFT {
				fnew -= a.gapextend
			} else {
				fnew -= a.gapopen
			}
		}
		if match > fnew && match > 0 {
			a.matrix[i][0] = match
			a.trace[i][0] = ALIGN_DIAG
		} else if fnew > 0 {
			a.matrix[i][0] = fnew
			a.trace[i][0] = ALIGN_LEFT
		} else {
			a.matrix[i][0] = 0.0
			a.trace[i][0] = ALIGN_DIAG // TO REVIEW
		}
	}

	// First column
	for j = 0; j <= l2; j++ {
		c1 = a.seq1.SequenceChar(0)
		c2 = a.seq2.SequenceChar(j)
		match = a.matchScore(c1, c2, indexseq1[0], indexseq2[j])
		fnew = 0.0
		if j > 0 {
			fnew = a.matrix[0][j-1]
			if a.trace[0][j-1] == ALIGN_UP {
				fnew -= a.gapextend
			} else {
				fnew -= a.gapopen
			}
		}
		if match > fnew && result > 0 {
			a.matrix[0][j] = match
			a.trace[0][j] = ALIGN_DIAG
		} else if fnew > 0 {
			a.matrix[0][j] = fnew
			a.trace[0][j] = ALIGN_UP
		} else {
			a.matrix[0][j] = 0.0
			a.trace[0][j] = ALIGN_DIAG // TO REVIEW
		}

		if j > 0 {
			maxa[j] = a.matrix[0][j]
			if a.trace[0][j-1] == ALIGN_UP {
				maxa[j] -= gapextend
			} else {
				maxa[j] -= gapopen
			}
		} else {
			maxa[j] = a.matrix[0][j] - a.gapopen
		}
	}

	for i = 1; i <= len(a.seq1.SequenceChar()); i++ {
		c1 = a.seq1.SequenceChar()[i-1]
		for j = 1; j <= len(a.seq2.SequenceChar()); j++ {
			c2 = a.seq2.SequenceChar()[j-1]

			maxscore = 0.0
			trace = ALIGN_STOP

			// left score
			leftscore = 0.0
			leftscore = a.matrix[i-1][j] + a.gap
			if leftscore > maxscore {
				trace = ALIGN_LEFT
				maxscore = leftscore
			}

			// up score
			upscore = 0.0
			upscore = a.matrix[i][j-1] + a.gap
			if upscore > maxscore {
				trace = ALIGN_UP
				maxscore = upscore
			}

			// diag score
			diagscore = a.matrix[i-1][j-1]
			diagscore += a.matchScore(c1, c2, indexseq1[i], indexseq2[j])

			if diagscore > maxscore {
				trace = ALIGN_DIAG
				maxscore = diagscore
			}

			a.matrix[i][j] = maxscore
			a.trace[i][j] = trace
			if maxscore > a.maxscore {
				a.maxi = i
				a.maxj = j
				a.maxscore = maxscore
			}
		}
	}
	return
}

func (a *pwaligner) fillMatrix_ATG() (err error) {
	var i, j int
	var revi, revj int
	var c1, c2 rune
	var maxscore float64
	var trace int
	var leftscore, upscore, diagscore float64
	var indexseq1, indexseq2 []int // convert characters to subst matrix positions
	a.initMatrix(a.seq1.Length(), a.seq2.Length())

	if a.submatrix != nil {
		if indexseq1, err = a.seqToindices(a.seq1); err != nil {
			return
		}
		if indexseq2, err = a.seqToindices(a.seq2); err != nil {
			return
		}
	}

	// In ATG algo, we look at sequences in reverse order
	for revi = len(a.seq1.SequenceChar()); revi > 0; revi-- {
		c1 = a.seq1.SequenceChar()[revi-1]
		i = len(a.seq1.SequenceChar()) - revi + 1
		for revj = len(a.seq2.SequenceChar()); revj > 0; revj-- {
			c2 = a.seq2.SequenceChar()[revj-1]
			j = len(a.seq2.SequenceChar()) - revj + 1
			maxscore = 0.0
			trace = ALIGN_STOP

			// left score
			leftscore = 0.0
			leftscore = a.matrix[i-1][j] + a.gap
			if leftscore > maxscore {
				trace = ALIGN_LEFT
				maxscore = leftscore
			}

			// up score
			upscore = 0.0
			upscore = a.matrix[i][j-1] + a.gap
			if upscore > maxscore {
				trace = ALIGN_UP
				maxscore = upscore
			}

			// diag score
			diagscore = a.matrix[i-1][j-1]
			if a.submatrix != nil {
				diagscore += a.submatrix[indexseq1[revi-1]][indexseq2[revj-1]]
			} else {
				if c1 != c2 {
					diagscore += a.mismatch
				} else {
					diagscore += a.match
				}
			}
			if diagscore > maxscore {
				trace = ALIGN_DIAG
				maxscore = diagscore
			}

			a.matrix[i][j] = maxscore
			a.trace[i][j] = trace

			// We update the max score only for the last column of the matrix
			// i.e. revi=0
			if maxscore > a.maxscore && revi == 1 {
				a.maxi = i
				a.maxj = j
				a.maxscore = maxscore
			}
		}
	}
	return
}

// Indices of alignment end
func (a *pwaligner) AlignEnds() (int, int) {
	return a.end1, a.end2
}

// Indices of alignment end
func (a *pwaligner) AlignStarts() (int, int) {
	return a.start1, a.start2
}

func (a *pwaligner) backTrack() {
	switch a.algo {
	case ALIGN_ALGO_ATG:
		a.backTrack_ATG()
	default:
		a.backTrack_SW()
	}
}

func (a *pwaligner) backTrack_SW() {
	var i, j int
	var seq1, seq2 []rune

	a.end1 = a.maxi - 1
	a.end2 = a.maxj - 1

	i = a.maxi
	j = a.maxj

	seq1 = make([]rune, 0, 20)
	seq2 = make([]rune, 0, 20)

	for a.trace[i][j] != ALIGN_STOP {
		switch a.trace[i][j] {
		case ALIGN_UP:
			seq1 = append(seq1, '-')
			seq2 = append(seq2, a.seq2.CharAt(j-1))
			j--

		case ALIGN_DIAG:
			seq1 = append(seq1, a.seq1.CharAt(i-1))
			seq2 = append(seq2, a.seq2.CharAt(j-1))
			i--
			j--
		case ALIGN_LEFT:
			seq1 = append(seq1, a.seq1.CharAt(i-1))
			seq2 = append(seq2, '-')
			i--
		}
	}

	a.start1 = i - 1
	a.start2 = j - 1

	Reverse(seq1)
	Reverse(seq2)

	a.seq1ali = seq1
	a.seq2ali = seq2
}

// In ATG backtrack we take the maximum found from the ATG of seq1,
// i.e its start, i.e. the last column of the matrix
func (a *pwaligner) backTrack_ATG() {
	var i, j int
	var revi, revj int
	var len1, len2 int
	var seq1, alistr, seq2 []rune

	len1 = len(a.seq1.SequenceChar())
	len2 = len(a.seq2.SequenceChar())

	a.start1 = len1 - a.maxi
	a.start2 = len2 - a.maxj

	revi = a.maxi
	revj = a.maxj

	i = len1 - revi
	j = len2 - revj

	seq1 = make([]rune, 0, 20)
	seq2 = make([]rune, 0, 20)
	alistr = make([]rune, 0, 20)

	for a.trace[revi][revj] != ALIGN_STOP {
		switch a.trace[revi][revj] {
		case ALIGN_UP:
			seq1 = append(seq1, '-')
			seq2 = append(seq2, a.seq2.CharAt(j))
			alistr = append(alistr, ' ')
			revj--
			j++

		case ALIGN_DIAG:
			seq1 = append(seq1, a.seq1.CharAt(i))
			seq2 = append(seq2, a.seq2.CharAt(j))
			if a.seq1.CharAt(i) == a.seq2.CharAt(j) {
				alistr = append(alistr, '|')
			} else {
				alistr = append(alistr, ' ')
			}
			revi--
			revj--
			i++
			j++
		case ALIGN_LEFT:
			seq1 = append(seq1, a.seq1.CharAt(i))
			seq2 = append(seq2, '-')
			alistr = append(alistr, ' ')
			revi--
			i++
		}
	}

	a.end1 = i
	a.end2 = j

	a.seq1ali = seq1
	a.seq2ali = seq2
	a.alistr = alistr
}

func (a *pwaligner) AlignmentStr() string {
	alistr := make([]rune, 0, len(a.seq1ali)+len(a.alistr)+len(a.seq2ali)+3)
	alistr = append(alistr, a.seq1ali...)
	alistr = append(alistr, '\n')
	alistr = append(alistr, a.alistr...)
	alistr = append(alistr, '\n')
	alistr = append(alistr, a.seq2ali...)
	alistr = append(alistr, '\n')
	return string(alistr)
}

func (a *pwaligner) Alignment() (align Alignment, err error) {
	err = a.fillMatrix()
	a.backTrack()
	align = NewAlign(UNKNOWN)
	align.AddSequenceChar(a.seq1.Name(), a.seq1ali, "")
	align.AddSequenceChar(a.seq2.Name(), a.seq2ali, "")
	align.AutoAlphabet()
	return
}

// Converts the aa/nt sequences into chartopos indices
// of the substitution matrix
func (a *pwaligner) seqToindices(s Sequence) (indices []int, err error) {
	var i int
	var ok bool

	indices = make([]int, s.Length())
	for i = 1; i < len(s.SequenceChar()); i++ {
		indices[i], ok = a.chartopos[s.CharAt(i)]
		if !ok {
			err = fmt.Errorf("Character not part of alphabet : %c", s.CharAt(i))
			return
		}
	}
	return
}

func (a *pwaligner) matchScore(c1, c2 rune, i1, i2 int) (score float64) {
	if a.submatrix != nil {
		score = a.submatrix[i1][i2]
	} else {
		if c1 != c2 {
			score = a.mismatch
		} else {
			score = a.match
		}
	}
	return
}
