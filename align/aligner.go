package align

type PairwiseAligner interface {
	AlignEnds()
	AlignStarts()
	SetGapStartScore(gap float64)
	SetGapElongScore(gap float64)
	SetMismatchScpre(mismatch float64)
	SetMatchScore(mismatch float64)
	Alignment() Alignment
}

const (
	ALIGN_UP = iota
	ALIGN_LEFT
	ALIGN_DIAG
	ALIGN_STOP
)

type swaligner struct {
	matrix [][]float64 /* trace matrix */
	trace  [][]int

	maxscore         float64 // Maximum score of the matrix
	maxi, maxj       int     // Indices of the maximum score
	start1, start2   int
	end1, end2       int
	seq1, seq2       Sequence
	seq1ali, seq2ali []rune // aligned sequences
	gapstart         float64
	gapelong         float64
	match            float64
	mismatch         float64
}

func NewSwAligner(seq1, seq2 Sequence) *swaligner {
	return &swaligner{
		matrix:   nil,
		trace:    nil,
		maxscore: 0.0,
		maxi:     0,
		maxj:     0,
		start1:   0,
		start2:   0,
		end1:     0,
		end2:     0,
		seq1:     seq1,
		seq2:     seq2,
		gapstart: -1.0,
		gapelong: -0.5,
		match:    1.0,
		mismatch: -1.0,
	}
}

func (a *swaligner) initMatrix(l1, l2 int) {
	var i int

	a.matrix = make([][]float64, l1)
	a.trace = make([][]int, l1)
	for i, _ = range a.matrix {
		a.matrix[i] = make([]float64, l2)
		a.trace[i] = make([]int, l2)
	}
}

func (a *swaligner) SetGapOpenScore(gap float64) {
	a.gapstart = gap
}
func (a *swaligner) SetGapElongScore(gap float64) {
	a.gapelong = gap
}

func (a *swaligner) SetMismatchScore(mismatch float64) {
	a.mismatch = mismatch
}

func (a *swaligner) SetMatchScore(match float64) {
	a.match = match
}

func (a *swaligner) fillMatrix() {
	var i, j int
	var c1, c2 rune
	var maxscore float64
	var trace int
	var leftscore, upscore, diagscore float64

	a.initMatrix(a.seq1.Length(), a.seq2.Length())

	for i, c1 = range a.seq1.SequenceChar() {
		for j, c2 = range a.seq2.SequenceChar() {
			maxscore = 0.0
			trace = ALIGN_STOP

			// left score
			leftscore = 0.0
			if i > 0 {
				if a.trace[i-1][j] == ALIGN_LEFT {
					leftscore = a.matrix[i-1][j] + a.gapelong
				} else {
					leftscore = a.matrix[i-1][j] + a.gapstart
				}
			}
			if leftscore > maxscore {
				trace = ALIGN_LEFT
				maxscore = leftscore
			}

			// up score
			upscore = 0.0
			if j > 0 {
				if a.trace[i][j] == ALIGN_UP {
					upscore = a.matrix[i][j-1] + a.gapelong
				} else {
					upscore = a.matrix[i][j-1] + a.gapstart
				}
			}
			if upscore > maxscore {
				trace = ALIGN_UP
				maxscore = upscore
			}

			// diag score
			diagscore = 0.0
			if i > 0 && j > 0 {
				diagscore = a.matrix[i-1][j-1]
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
			if maxscore > a.maxscore {
				a.maxi = i
				a.maxj = j
				a.maxscore = maxscore
			}
		}
	}
}

// Indices of alignment end
func (a *swaligner) AlignEnds() (int, int) {
	return a.end1, a.end2
}

// Indices of alignment end
func (a *swaligner) AlignStarts() (int, int) {
	return a.start2, a.start2
}

func (a *swaligner) backTrack() {
	var i, j int
	var seq1, seq2 []rune

	a.end1 = a.maxi
	a.end2 = a.maxj

	i = a.end1
	j = a.end2

	seq1 = make([]rune, 0, 20)
	seq2 = make([]rune, 0, 20)

	for a.matrix[i][j] != 0.0 {
		switch a.trace[i][j] {
		case ALIGN_UP:
			seq1 = append(seq1, '-')
			seq2 = append(seq2, a.seq2.CharAt(j))
			j--
		case ALIGN_DIAG:
			seq1 = append(seq1, a.seq1.CharAt(i))
			seq2 = append(seq2, a.seq2.CharAt(j))
			i--
			j--
		case ALIGN_LEFT:
			seq1 = append(seq1, a.seq1.CharAt(i))
			seq2 = append(seq2, '-')
			i--
		}
	}

	a.start1 = i
	a.start2 = j

	Reverse(seq1)
	Reverse(seq2)

	a.seq1ali = seq1
	a.seq2ali = seq2
}

func (a *swaligner) Alignment() Alignment {
	a.fillMatrix()
	a.backTrack()
	align := NewAlign(UNKNOWN)
	align.AddSequenceChar(a.seq1.Name(), a.seq1ali, "")
	align.AddSequenceChar(a.seq2.Name(), a.seq2ali, "")
	align.AutoAlphabet()
	return align
}
