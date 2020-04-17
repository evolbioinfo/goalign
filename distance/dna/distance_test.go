package dna

import (
	"math"
	"testing"

	"github.com/evolbioinfo/goalign/align"
)

func Test_countDiffs(t *testing.T) {
	selectedSites := []bool{true, true, true, true, true, true, true, true, true, true, true}
	weights := []float64{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}
	type args struct {
		seq1          []rune
		seq2          []rune
		selectedSites []bool
		weights       []float64
	}
	tests := []struct {
		name        string
		args        args
		wantNbdiffs float64
		wantTotal   float64
	}{
		{name: "t1", args: args{seq1: []rune("ACGTACGTNNR"), seq2: []rune("ACGTACGTNNR"), selectedSites: selectedSites, weights: weights}, wantNbdiffs: 0.0, wantTotal: 11.0},
		{name: "t2", args: args{seq1: []rune("ACGTACGTNNR"), seq2: []rune("ACCTACGTNNR"), selectedSites: selectedSites, weights: weights}, wantNbdiffs: 1.0, wantTotal: 11.0},
		{name: "t3", args: args{seq1: []rune("ACGT-CGTNNR"), seq2: []rune("ACCTACGTNNR"), selectedSites: selectedSites, weights: weights}, wantNbdiffs: 1.0, wantTotal: 10.0},
		{name: "t4", args: args{seq1: []rune("ACGTACGTNNR"), seq2: []rune("ACGTACGTNNS"), selectedSites: selectedSites, weights: weights}, wantNbdiffs: 2.0 / 3.0, wantTotal: 11.0},
		{name: "t5", args: args{seq1: []rune("ACGTACGTNN-"), seq2: []rune("ACGTACGTNNR"), selectedSites: selectedSites, weights: weights}, wantNbdiffs: 0.0, wantTotal: 10.0},
		{name: "t6", args: args{seq1: []rune("ACGTACGTNNR"), seq2: []rune("ACGTACGTNGR"), selectedSites: selectedSites, weights: weights}, wantNbdiffs: 3.0 / 4.0, wantTotal: 11.0},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			var s1, s2 []uint8
			var err error

			s1 = make([]uint8, len(tt.args.seq1))
			for l, r := range tt.args.seq1 {
				if s1[l], err = align.Nt2IndexIUPAC(r); err != nil {
					t.Error(err)
				}
			}
			s2 = make([]uint8, len(tt.args.seq2))
			for l, r := range tt.args.seq2 {
				if s2[l], err = align.Nt2IndexIUPAC(r); err != nil {
					t.Error(err)
				}
			}

			gotNbdiffs, gotTotal := countDiffs(s1, s2, tt.args.selectedSites, tt.args.weights)
			if math.Abs(gotNbdiffs-tt.wantNbdiffs) > 0.000000000000001 {
				t.Errorf("countDiffs() gotNbdiffs = %v, want %v", gotNbdiffs, tt.wantNbdiffs)
			}
			if math.Abs(gotTotal-tt.wantTotal) > 0.000000000000001 {
				t.Errorf("countDiffs() gotTotal = %v, want %v", gotTotal, tt.wantTotal)
			}
		})
	}
}

func Test_countDiffsWithGaps(t *testing.T) {
	selectedSites := []bool{true, true, true, true, true, true, true, true, true, true, true}
	weights := []float64{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}

	type args struct {
		seq1          []rune
		seq2          []rune
		selectedSites []bool
		weights       []float64
	}
	tests := []struct {
		name        string
		args        args
		wantNbdiffs float64
		wantTotal   float64
	}{
		{name: "t1", args: args{seq1: []rune("ACGTACGTNNR"), seq2: []rune("ACGTACGTNNR"), selectedSites: selectedSites, weights: weights}, wantNbdiffs: 0.0, wantTotal: 11.0},
		{name: "t2", args: args{seq1: []rune("ACGTACGTNNR"), seq2: []rune("ACCTACGTNNR"), selectedSites: selectedSites, weights: weights}, wantNbdiffs: 1.0, wantTotal: 11.0},
		{name: "t3", args: args{seq1: []rune("ACGT-CGTNNR"), seq2: []rune("ACCTACGTNNR"), selectedSites: selectedSites, weights: weights}, wantNbdiffs: 2.0, wantTotal: 11.0},
		{name: "t4", args: args{seq1: []rune("ACGTACGTNNR"), seq2: []rune("ACGTACGTNNS"), selectedSites: selectedSites, weights: weights}, wantNbdiffs: 2.0 / 3.0, wantTotal: 11.0},
		{name: "t5", args: args{seq1: []rune("ACGTACGTNN-"), seq2: []rune("ACGTACGTNNR"), selectedSites: selectedSites, weights: weights}, wantNbdiffs: 1.0, wantTotal: 11.0},
		{name: "t6", args: args{seq1: []rune("ACGTACGTNNR"), seq2: []rune("ACGTACGTNGR"), selectedSites: selectedSites, weights: weights}, wantNbdiffs: 3.0 / 4.0, wantTotal: 11.0},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			var s1, s2 []uint8
			var err error

			s1 = make([]uint8, len(tt.args.seq1))
			for l, r := range tt.args.seq1 {
				if s1[l], err = align.Nt2IndexIUPAC(r); err != nil {
					t.Error(err)
				}
			}
			s2 = make([]uint8, len(tt.args.seq2))
			for l, r := range tt.args.seq2 {
				if s2[l], err = align.Nt2IndexIUPAC(r); err != nil {
					t.Error(err)
				}
			}

			gotNbdiffs, gotTotal := countDiffsWithGaps(s1, s2, tt.args.selectedSites, tt.args.weights)
			if math.Abs(gotNbdiffs-tt.wantNbdiffs) > 0.000000000000001 {
				t.Errorf("countDiffs() gotNbdiffs = %v, want %v", gotNbdiffs, tt.wantNbdiffs)
			}
			if math.Abs(gotTotal-tt.wantTotal) > 0.000000000000001 {
				t.Errorf("countDiffs() gotTotal = %v, want %v", gotTotal, tt.wantTotal)
			}
		})
	}
}

func Test_countDiffsWithInternalGaps(t *testing.T) {
	selectedSites := []bool{true, true, true, true, true, true, true, true, true, true, true}
	weights := []float64{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}

	type args struct {
		seq1          []rune
		seq2          []rune
		selectedSites []bool
		weights       []float64
	}
	tests := []struct {
		name        string
		args        args
		wantNbdiffs float64
		wantTotal   float64
	}{
		{name: "t1", args: args{seq1: []rune("ACGTACGTNNR"), seq2: []rune("ACGTACGTNNR"), selectedSites: selectedSites, weights: weights}, wantNbdiffs: 0.0, wantTotal: 11.0},
		{name: "t2", args: args{seq1: []rune("ACGTACGTNNR"), seq2: []rune("ACCTACGTNNR"), selectedSites: selectedSites, weights: weights}, wantNbdiffs: 1.0, wantTotal: 11.0},
		{name: "t3", args: args{seq1: []rune("ACGT-CGTNNR"), seq2: []rune("ACCTACGTNNR"), selectedSites: selectedSites, weights: weights}, wantNbdiffs: 2.0, wantTotal: 11.0},
		{name: "t4", args: args{seq1: []rune("ACGTACGTNNR"), seq2: []rune("ACGTACGTNNS"), selectedSites: selectedSites, weights: weights}, wantNbdiffs: 2.0 / 3.0, wantTotal: 11.0},
		{name: "t5", args: args{seq1: []rune("-CGTACGTNN-"), seq2: []rune("ACGTACGTNNR"), selectedSites: selectedSites, weights: weights}, wantNbdiffs: 0.0, wantTotal: 9.0},
		{name: "t6", args: args{seq1: []rune("ACGTACGTNNR"), seq2: []rune("ACGTACGTNGR"), selectedSites: selectedSites, weights: weights}, wantNbdiffs: 3.0 / 4.0, wantTotal: 11.0},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			var s1, s2 []uint8
			var err error

			s1 = make([]uint8, len(tt.args.seq1))
			for l, r := range tt.args.seq1 {
				if s1[l], err = align.Nt2IndexIUPAC(r); err != nil {
					t.Error(err)
				}
			}
			s2 = make([]uint8, len(tt.args.seq2))
			for l, r := range tt.args.seq2 {
				if s2[l], err = align.Nt2IndexIUPAC(r); err != nil {
					t.Error(err)
				}
			}

			gotNbdiffs, gotTotal := countDiffsWithInternalGaps(s1, s2, tt.args.selectedSites, tt.args.weights)
			if math.Abs(gotNbdiffs-tt.wantNbdiffs) > 0.000000000000001 {
				t.Errorf("countDiffs() gotNbdiffs = %v, want %v", gotNbdiffs, tt.wantNbdiffs)
			}
			if math.Abs(gotTotal-tt.wantTotal) > 0.000000000000001 {
				t.Errorf("countDiffs() gotTotal = %v, want %v", gotTotal, tt.wantTotal)
			}
		})
	}
}
