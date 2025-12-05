package dna

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
)

// Compute the distance for each windows between the reference sequence and all the
// other sequences. Can be used to draw a simplot
func SimPlotDistances(al align.Alignment, refseq string, distmodel string, windowsize, windowstep int) (sp []struct {
	WindowStart, WindowEnd int
	CompSeq                string
	Distance               float64
}, err error) {
	var subal align.Alignment
	var model DistModel

	// We create the dist windows array
	// with size (nseq-1)*nwindows
	sp = make([]struct {
		WindowStart int
		WindowEnd   int
		CompSeq     string
		Distance    float64
	}, 0, (al.NbSequences()-1)*(al.Length()/windowsize))

	if al.Alphabet() != align.NUCLEOTIDS {
		err = fmt.Errorf("alignment must be nucleotidic")
		return
	}

	// Extract reference sequence id (sequence to compare to all)
	refid := al.GetSequenceIdByName(refseq)
	if refid < 0 {
		err = fmt.Errorf("reference sequence %s does not exist in the alignment", refseq)
		return
	}

	// Initialize model
	if model, err = Model(distmodel, false); err != nil {
		io.LogError(err)
		return
	}

	// Iterate over the windows
	for start := 0; start+windowsize <= al.Length(); start += windowstep {

		if subal, err = al.SubAlign(start, windowsize); err != nil {
			return
		}

		var range1min, range1max, range2min, range2max int = refid, refid, -1, -1

		// Output matrix for this sub alignment
		var distMatrix [][]float64
		distMatrix, err = DistMatrix(subal, nil, model, range1min, range1max, range2min, range2max, false, 0.0, 1)
		if err != nil {
			io.LogError(err)
			return
		}

		for i := 0; i < al.NbSequences(); i++ {
			if i != refid {
				compname, _ := al.GetSequenceNameById(i)
				window := struct {
					WindowStart, WindowEnd int
					CompSeq                string
					Distance               float64
				}{start, start + windowsize, compname, distMatrix[refid][i]}
				sp = append(sp, window)
			}
		}
	}
	return
}
