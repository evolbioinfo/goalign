package align

import (
	"math"
	"reflect"
	"testing"
)

func TestGenAllPossibleCodons(t *testing.T) {
	type args struct {
		nt1 rune
		nt2 rune
		nt3 rune
	}
	tests := []struct {
		name       string
		args       args
		wantCodons []string
	}{
		{name: "AGN", args: args{nt1: 'A', nt2: 'G', nt3: 'N'}, wantCodons: []string{"AGA", "AGC", "AGG", "AGT"}},
		{name: "AGZ", args: args{nt1: 'A', nt2: 'G', nt3: 'Z'}, wantCodons: []string{}},
		{name: "AGT", args: args{nt1: 'A', nt2: 'G', nt3: 'T'}, wantCodons: []string{"AGT"}},
		{name: "ANT", args: args{nt1: 'A', nt2: 'N', nt3: 'T'}, wantCodons: []string{"AAT", "ACT", "AGT", "ATT"}},
		{name: "ANN", args: args{nt1: 'A', nt2: 'N', nt3: 'N'}, wantCodons: []string{"AAA", "ACA", "AGA", "ATA", "AAC", "ACC", "AGC", "ATC", "AAG", "ACG", "AGG", "ATG", "AAT", "ACT", "AGT", "ATT"}},
		{name: "AGT", args: args{nt1: 'A', nt2: 'G', nt3: 'T'}, wantCodons: []string{"AGT"}},
		{name: "AYG", args: args{nt1: 'A', nt2: 'Y', nt3: 'G'}, wantCodons: []string{"ACG", "ATG"}},
		{name: "RGA", args: args{nt1: 'R', nt2: 'G', nt3: 'A'}, wantCodons: []string{"AGA", "GGA"}},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if gotCodons := GenAllPossibleCodons(tt.args.nt1, tt.args.nt2, tt.args.nt3); !reflect.DeepEqual(gotCodons, tt.wantCodons) {
				t.Errorf("GenAllPossibleCodons() = %v, want %v", gotCodons, tt.wantCodons)
			}
		})
	}
}

func Test_seq_Translate(t *testing.T) {
	type fields struct {
		name     string
		sequence []rune
		comment  string
	}
	type args struct {
		phase       int
		geneticcode int
	}
	tests := []struct {
		name    string
		fields  fields
		args    args
		wantTr  Sequence
		wantErr bool
	}{
		{name: "Seq1", fields: fields{name: "seq1", sequence: []rune{'G', 'A', 'Y', 'A', 'A', 'R', 'U', 'A', 'Y', 'C', 'A', 'Y', 'R', 'A', 'Y', 'U', 'A', 'G'}}, args: args{phase: 0, geneticcode: 0}, wantTr: &seq{name: "seq1", sequence: []rune{'D', 'K', 'Y', 'H', 'X', '*'}}, wantErr: false},
		{name: "Seq1", fields: fields{name: "seq1", sequence: []rune{'G', 'A', 'Y', 'A', 'A', 'R', 'U', 'A', 'Y', 'C', 'A', 'Y', 'A', 'A', 'Y', 'U', 'A', 'G'}}, args: args{phase: 0, geneticcode: 0}, wantTr: &seq{name: "seq1", sequence: []rune{'D', 'K', 'Y', 'H', 'N', '*'}}, wantErr: false},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			s := &seq{
				name:     tt.fields.name,
				sequence: tt.fields.sequence,
				comment:  tt.fields.comment,
			}
			gotTr, err := s.Translate(tt.args.phase, tt.args.geneticcode)
			if (err != nil) != tt.wantErr {
				t.Errorf("seq.Translate() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(gotTr, tt.wantTr) {
				t.Errorf("seq.Translate() = %v, want %v", gotTr, tt.wantTr)
			}
		})
	}
}

func Test_seq_NumGapsFromEnd(t *testing.T) {
	type fields struct {
		name     string
		sequence []rune
		comment  string
	}
	tests := []struct {
		name        string
		fields      fields
		wantNumgaps int
	}{
		{name: "num 0", fields: fields{name: "s1", sequence: []rune{'-', '-', '-', 'A', '-', 'C'}, comment: ""}, wantNumgaps: 0},
		{name: "num 1", fields: fields{name: "s1", sequence: []rune{'G', '-', '-', '-', 'A', '-', '-'}, comment: ""}, wantNumgaps: 2},
		{name: "num 2", fields: fields{name: "s1", sequence: []rune{'-', '-', '-', '-', '-', '-'}, comment: ""}, wantNumgaps: 6},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			s := &seq{
				name:     tt.fields.name,
				sequence: tt.fields.sequence,
				comment:  tt.fields.comment,
			}
			if gotNumgaps := s.NumGapsFromEnd(); gotNumgaps != tt.wantNumgaps {
				t.Errorf("seq.NumGapsFromEnd() = %v, want %v", gotNumgaps, tt.wantNumgaps)
			}
		})
	}
}

func Test_seq_NumGapsFromStart(t *testing.T) {
	type fields struct {
		name     string
		sequence []rune
		comment  string
	}
	tests := []struct {
		name        string
		fields      fields
		wantNumgaps int
	}{
		{name: "num 0", fields: fields{name: "s1", sequence: []rune{'-', '-', '-', 'A', '-', 'C'}, comment: ""}, wantNumgaps: 3},
		{name: "num 1", fields: fields{name: "s1", sequence: []rune{'G', '-', '-', '-', 'A', '-', 'C'}, comment: ""}, wantNumgaps: 0},
		{name: "num 2", fields: fields{name: "s1", sequence: []rune{'-', '-', '-', '-', '-', '-'}, comment: ""}, wantNumgaps: 6},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			s := &seq{
				name:     tt.fields.name,
				sequence: tt.fields.sequence,
				comment:  tt.fields.comment,
			}
			if gotNumgaps := s.NumGapsFromStart(); gotNumgaps != tt.wantNumgaps {
				t.Errorf("seq.NumGapsFromStart() = %v, want %v", gotNumgaps, tt.wantNumgaps)
			}
		})
	}
}

func Test_seq_NumGaps(t *testing.T) {
	type fields struct {
		name     string
		sequence []rune
		comment  string
	}
	tests := []struct {
		name        string
		fields      fields
		wantNumgaps int
	}{
		{name: "num 0", fields: fields{name: "s1", sequence: []rune{'-', '-', '-', 'A', '-', 'C'}, comment: ""}, wantNumgaps: 4},
		{name: "num 1", fields: fields{name: "s1", sequence: []rune{'G', '-', '-', '-', 'A', '-', 'C'}, comment: ""}, wantNumgaps: 4},
		{name: "num 2", fields: fields{name: "s1", sequence: []rune{'-', '-', '-', '-', '-', '-'}, comment: ""}, wantNumgaps: 6},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			s := &seq{
				name:     tt.fields.name,
				sequence: tt.fields.sequence,
				comment:  tt.fields.comment,
			}
			if gotNumgaps := s.NumGaps(); gotNumgaps != tt.wantNumgaps {
				t.Errorf("seq.NumGaps() = %v, want %v", gotNumgaps, tt.wantNumgaps)
			}
		})
	}
}

func TestEqualOrCompatible(t *testing.T) {
	type args struct {
		nt1 rune
		nt2 rune
	}
	tests := []struct {
		name    string
		args    args
		wantOk  bool
		wantErr bool
	}{
		{name: "t1", args: args{nt1: 'S', nt2: 'Y'}, wantOk: true, wantErr: false},
		{name: "t2", args: args{nt1: 'R', nt2: 'Y'}, wantOk: false, wantErr: false},
		{name: "t3", args: args{nt1: 'W', nt2: 'K'}, wantOk: true, wantErr: false},
		{name: "t4", args: args{nt1: 'N', nt2: 'A'}, wantOk: true, wantErr: false},
		{name: "t5", args: args{nt1: 'A', nt2: 'C'}, wantOk: false, wantErr: false},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotOk, err := EqualOrCompatible(tt.args.nt1, tt.args.nt2)
			if (err != nil) != tt.wantErr {
				t.Errorf("EqualOrCompatible() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if gotOk != tt.wantOk {
				t.Errorf("EqualOrCompatible() = %v, want %v", gotOk, tt.wantOk)
			}
		})
	}
}

func TestNtIUPACDifference(t *testing.T) {
	type args struct {
		nt1 rune
		nt2 rune
	}
	tests := []struct {
		name     string
		args     args
		wantDiff float64
		wantErr  bool
	}{
		{name: "t1", args: args{nt1: 'S', nt2: 'Y'}, wantDiff: 1.0 - 1.0/3.0, wantErr: false},
		{name: "t2", args: args{nt1: 'R', nt2: 'Y'}, wantDiff: 1.0, wantErr: false},
		{name: "t3", args: args{nt1: 'W', nt2: 'K'}, wantDiff: 1.0 - 1.0/3.0, wantErr: false},
		{name: "t4", args: args{nt1: 'N', nt2: 'A'}, wantDiff: 3.0 / 4.0, wantErr: false},
		{name: "t5", args: args{nt1: 'A', nt2: 'C'}, wantDiff: 1.0, wantErr: false},
		{name: "t6", args: args{nt1: 'C', nt2: 'C'}, wantDiff: 0.0, wantErr: false},
		{name: "t7", args: args{nt1: '-', nt2: 'C'}, wantDiff: 1.0, wantErr: false},
		{name: "t8", args: args{nt1: 'C', nt2: '-'}, wantDiff: 1.0, wantErr: false},
		{name: "t9", args: args{nt1: '-', nt2: '-'}, wantDiff: 0.0, wantErr: false},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			gotDiff, err := NtIUPACDifference(tt.args.nt1, tt.args.nt2)
			if (err != nil) != tt.wantErr {
				t.Errorf("NtIUPACDifference() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if math.Abs(gotDiff-tt.wantDiff) > 0.000000000000001 {
				t.Errorf("NtIUPACDifference() = %v, want %v", gotDiff, tt.wantDiff)
			}
		})
	}
}
