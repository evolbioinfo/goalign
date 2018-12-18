package align

import (
	"fmt"
)

type PairwiseAligner interface {
	AlignEnds() (int, int)
	AlignStarts() (int, int)
	Seq1Ali() []rune
	Seq2Ali() []rune
	SetGapOpenScore(open float64)
	SetGapExtendScore(extend float64)
	SetScore(match, mismatch float64)
	MaxScore() float64 // Maximum score of the alignment
	NbMatches() int    // Number of matches
	NbMisMatches() int // Number of mismatches
	NbGaps() int       // Nuber of gaps
	Length() int       // Length of the alignment
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

type pwaligner struct {
	algo   int
	matrix [][]float64 // dynamix programming matrix
	trace  [][]int     // trace matrix
	maxa   []float64   // keep track of best gap opened

	maxscore         float64 // Maximum score of the matrix
	nbmatches        int     // number of matches
	nbmismatches     int     // number of mismatches
	nbgaps           int     // number of gaps
	length           int     // alignment length
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

	if (a1 == NUCLEOTIDS || a1 == BOTH) && (a2 == NUCLEOTIDS || a2 == BOTH) {
		mat = dnafull_subst_matrix
		chartopos = dna_to_matrix_pos
	} else if (a1 == AMINOACIDS || a1 == BOTH) && (a2 == AMINOACIDS || a2 == BOTH) {
		mat = blosum62_subst_matrix
		chartopos = prot_to_matrix_pos
	}

	return &pwaligner{
		algo:      algo,
		matrix:    nil,
		trace:     nil,
		maxa:      nil,
		maxscore:  .0,
		maxi:      0,
		maxj:      0,
		start1:    0,
		start2:    0,
		end1:      0,
		end2:      0,
		seq1:      seq1.Clone(),
		seq2:      seq2.Clone(),
		gapopen:   -10.0,
		gapextend: -0.5,
		match:     1.0,
		mismatch:  -1.0,
		submatrix: mat,
		chartopos: chartopos,
	}
}

func (a *pwaligner) initMatrix(l1, l2 int) {
	var i int

	a.matrix = make([][]float64, l1)
	a.trace = make([][]int, l1)
	a.maxa = make([]float64, l2)
	for i, _ = range a.matrix {
		a.matrix[i] = make([]float64, l2)
		a.trace[i] = make([]int, l2)
	}
}

func (a *pwaligner) SetGapOpenScore(gap float64) {
	a.gapopen = gap
}

func (a *pwaligner) SetGapExtendScore(gap float64) {
	a.gapextend = gap
}

// Sets manually match and mismatch scores
// Substitution matrix is not used any more
func (a *pwaligner) SetScore(match, mismatch float64) {
	a.match = match
	a.mismatch = mismatch
	a.submatrix = nil
	//a.chartopos = nil
}

func (a *pwaligner) MaxScore() float64 {
	return a.maxscore
}

func (a *pwaligner) NbMatches() int {
	return a.nbmatches
}

func (a *pwaligner) NbMisMatches() int {
	return a.nbmismatches
}

func (a *pwaligner) NbGaps() int {
	return a.nbgaps
}

func (a *pwaligner) Length() int {
	return a.length
}

func (a *pwaligner) fillMatrix() (err error) {
	switch a.algo {
	case ALIGN_ALGO_ATG:
		// We want to backtrack from the atg of the first sequence
		// i.e. the max value of the last row
		a.seq1.Reverse()
		a.seq2.Reverse()
		err = a.fillMatrix_SW()
	default:
		err = a.fillMatrix_SW()
	}
	return
}

// Inspired from function embAlignPathCalcSW
// of EMBOSS suite
func (a *pwaligner) fillMatrix_SW() (err error) {
	var l1, l2 int
	var c1, c2 rune
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

	l1 = a.seq1.Length()
	l2 = a.seq2.Length()

	// We initialize first row and first column of the matrix
	var match, fnew float64

	// First row
	for j := 0; j < l2; j++ {
		c1 = a.seq1.CharAt(0)
		c2 = a.seq2.CharAt(j)
		match = a.matchScore(c1, c2, indexseq1[0], indexseq2[j])
		fnew = 0.0
		if j > 0 {
			fnew = a.matrix[0][j-1]
			if a.trace[0][j-1] == ALIGN_LEFT {
				fnew += a.gapextend
			} else {
				fnew += a.gapopen
			}
		}
		if match > fnew && match > .0 {
			a.matrix[0][j] = match
			a.trace[0][j] = ALIGN_DIAG
		} else if fnew > .0 {
			a.matrix[0][j] = fnew
			a.trace[0][j] = ALIGN_LEFT
		} else {
			a.matrix[0][j] = 0.0
			a.trace[0][j] = ALIGN_DIAG // TO REVIEW
		}

		if j > 0 {
			a.maxa[j] = a.matrix[0][j]
			if a.trace[0][j-1] == ALIGN_LEFT {
				a.maxa[j] += a.gapextend
			} else {
				a.maxa[j] += a.gapopen
			}
		} else {
			a.maxa[j] = a.matrix[0][j] + a.gapopen
		}
	}

	// First column
	for i := 0; i < l1; i++ {
		c1 = a.seq1.CharAt(i)
		c2 = a.seq2.CharAt(0)
		match = a.matchScore(c1, c2, indexseq1[i], indexseq2[0])

		fnew = 0.0
		if i > 0 {
			fnew = a.matrix[i-1][0]
			if a.trace[i-1][0] == ALIGN_UP {
				fnew += a.gapextend
			} else {
				fnew += a.gapopen
			}
		}
		if match > fnew && match > .0 {
			a.matrix[i][0] = match
			a.trace[i][0] = ALIGN_DIAG
		} else if fnew > 0 {
			a.matrix[i][0] = fnew
			a.trace[i][0] = ALIGN_UP
		} else {
			a.matrix[i][0] = 0.0
			a.trace[i][0] = ALIGN_DIAG // TO REVIEW
		}
	}

	// Each line
	for i := 1; i < l1; i++ {
		c1 = a.seq1.SequenceChar()[i]
		// Temp value for max left gap extensions of this line
		bx := a.matrix[i][0] + a.gapopen + a.gapextend
		// Each column
		for j := 1; j < l2; j++ {
			c2 = a.seq2.SequenceChar()[j]
			match = a.matchScore(c1, c2, indexseq1[i], indexseq2[j])

			// diag score
			mscore := a.matrix[i-1][j-1] + match
			a.trace[i][j] = ALIGN_DIAG
			a.matrix[i][j] = mscore

			// Search if best options with gaps
			// i direction
			// either gap opening or gap extention
			// without re looping over 0-i
			a.maxa[j] += a.gapextend
			fnew = a.matrix[i-1][j]
			fnew += a.gapopen
			if fnew > a.maxa[j] {
				a.maxa[j] = fnew
			}
			if a.maxa[j] > mscore {
				mscore = a.maxa[j]
				a.matrix[i][j] = mscore
				a.trace[i][j] = ALIGN_UP
			}

			// Search if best options with gaps
			// j direction
			// either gap opening or gap extention
			// without re looping over 0-j
			bx += a.gapextend
			fnew = a.matrix[i][j-1]
			fnew += a.gapopen
			if fnew > bx {
				bx = fnew
			}
			if bx > mscore {
				mscore = bx
				a.matrix[i][j] = mscore
				a.trace[i][j] = ALIGN_LEFT
			}
			if mscore > a.maxscore {
				a.maxscore = mscore
				a.maxi = i
				a.maxj = j
			}
			if a.matrix[i][j] < .0 {
				a.matrix[i][j] = .0
			}
		}
	}

	// for i, v := range a.matrix {
	// 	for j, v2 := range v {
	// 		fmt.Printf("\t%.2f", v2)
	// 		if a.matrix[i][j] == ALIGN_UP {
	// 			fmt.Printf("^")
	// 		} else if a.matrix[i][j] == ALIGN_LEFT {
	// 			fmt.Printf("<")
	// 		}
	// 	}
	// 	fmt.Println()
	// }

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
		// We want to backtrack from the atg of the first sequence
		// i.e. the max value of the last row
		// We set the max of last row
		a.maxi = 0
		a.maxj = 0
		a.maxscore = .0
		for j := 0; j < a.seq2.Length(); j++ {
			if a.matrix[a.seq1.Length()-1][j] > a.maxscore {
				a.maxscore = a.matrix[a.seq1.Length()-1][j]
				a.maxi = a.seq1.Length() - 1
				a.maxj = j
			}
		}
		a.backTrack_SW()
		// And we reverse the sequences back
		a.seq1.Reverse()
		a.seq2.Reverse()
		Reverse(a.alistr)
		Reverse(a.seq1ali)
		Reverse(a.seq2ali)
		a.start1, a.end1 = a.seq1.Length()-a.end1-1, a.seq1.Length()-a.start1-1
		a.start2, a.end2 = a.seq2.Length()-a.end2-1, a.seq2.Length()-a.start2-1
	default:
		a.backTrack_SW()
	}
}

func (a *pwaligner) backTrack_SW() {
	var i, j, ngaps int
	var seq1, seq2, alistr []rune
	var gapscore float64

	a.end1 = a.maxi
	a.end2 = a.maxj

	i = a.maxi
	j = a.maxj

	seq1 = make([]rune, 0, 20)
	seq2 = make([]rune, 0, 20)
	alistr = make([]rune, 0, 20)

	for i >= 0 && j >= 0 {
		switch a.trace[i][j] {
		case ALIGN_UP:
			ngaps = 0
			for {
				ngaps++
				gapscore = a.matrix[i-ngaps][j] + a.gapopen + float64(ngaps-1)*a.gapextend
				if gapscore == a.matrix[i][j] || i-ngaps == 0 {
					break
				}
			}
			for g := 0; g < ngaps; g++ {
				seq1 = append(seq1, a.seq1.CharAt(i))
				seq2 = append(seq2, '-')
				alistr = append(alistr, ' ')
				a.length++
				a.nbgaps++
				i--
			}
		case ALIGN_DIAG:
			seq1 = append(seq1, a.seq1.CharAt(i))
			seq2 = append(seq2, a.seq2.CharAt(j))
			a.length++
			if a.seq2.CharAt(j) == a.seq1.CharAt(i) {
				a.nbmatches++
				alistr = append(alistr, '|')
			} else {
				a.nbmismatches++
				alistr = append(alistr, '.')
			}
			i--
			j--
		case ALIGN_LEFT:
			ngaps = 0
			for {
				ngaps++
				gapscore = a.matrix[i][j-ngaps] + a.gapopen + float64(ngaps-1)*a.gapextend
				if gapscore == a.matrix[i][j] || j-ngaps == 0 {
					break
				}
			}
			for g := 0; g < ngaps; g++ {
				seq1 = append(seq1, '-')
				seq2 = append(seq2, a.seq2.CharAt(j))
				alistr = append(alistr, ' ')
				a.length++
				a.nbgaps++
				j--
			}
		}
		if i > 0 && j > 0 && a.matrix[i][j] <= .0 {
			break
		}
	}

	a.start1 = i + 1
	a.start2 = j + 1
	Reverse(seq1)
	Reverse(seq2)
	Reverse(alistr)

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

func (a *pwaligner) Seq1Ali() []rune {
	return a.seq1ali
}

func (a *pwaligner) Seq2Ali() []rune {
	return a.seq2ali
}

func (a *pwaligner) Alignment() (align Alignment, err error) {

	if err = a.fillMatrix(); err != nil {
		return
	}
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
	for i = 0; i < len(s.SequenceChar()); i++ {
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
