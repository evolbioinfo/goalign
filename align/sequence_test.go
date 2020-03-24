package align

import (
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
