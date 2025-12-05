package dna

import (
	"fmt"
	"strings"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"
)

// Compute the distance for each windows between the reference sequence and all the
// other sequences. Can be used to draw a simplot
// If group is true
// Then extract group name from the sequence names, using splitfield and splitsep
// Sequences are grouped and distance from refseq is computed as the average with each group
func SimPlotDistances(al align.Alignment, refseq string, distmodel string, windowsize, windowstep int, group bool, splitsep string, splitfield int) (sp []struct {
	WindowStart, WindowEnd int
	CompSeq                string
	Distance               float64
}, err error) {
	var subal align.Alignment
	var model DistModel
	// Groups of sequence ids
	var groups map[string][]int = make(map[string][]int)

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

	// We create groups
	if group {
		if splitfield < 0 {
			err = fmt.Errorf("groupe separator fiels must be positive")
			return
		}
		for i, s := range al.Sequences() {
			if s.Name() != refseq {
				cols := strings.Split(s.Name(), splitsep)
				if len(cols) < 2 && len(cols) <= splitfield {
					err = fmt.Errorf("split sequence name has less than 2 values or less than the desired field: %s", s.Name())
					return
				}
				gName := cols[splitfield]
				groups[gName] = append(groups[gName], i)
			}
		}
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

		if group {
			for gName := range groups {
				avg := 0.0
				for _, id := range groups[gName] {
					avg += distMatrix[refid][id]
				}
				avg /= float64(len(groups[gName]))
				window := struct {
					WindowStart, WindowEnd int
					CompSeq                string
					Distance               float64
				}{start, start + windowsize, gName, avg}
				sp = append(sp, window)
			}
		} else {
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
	}
	return
}
