package align

import (
	"reflect"
	"testing"
)

func Test_seqbag_UniqueCharacters(t *testing.T) {
	type fields struct {
		seqmap          map[string]*seq
		seqs            []*seq
		ignoreidentical int
		alphabet        int
	}
	tests := []struct {
		name      string
		fields    fields
		wantChars []uint8
	}{
		{name: "t1",
			fields: fields{seqmap: nil,
				seqs: []*seq{
					{sequence: []uint8("ACGTACGTACGT")},
					{sequence: []uint8("ACGTAC*TACGT")}}},
			wantChars: []uint8("*ACGT")},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			sb := &seqbag{
				seqmap:          tt.fields.seqmap,
				seqs:            tt.fields.seqs,
				ignoreidentical: tt.fields.ignoreidentical,
				alphabet:        tt.fields.alphabet,
			}
			if gotChars := sb.UniqueCharacters(); !reflect.DeepEqual(gotChars, tt.wantChars) {
				t.Errorf("seqbag.UniqueCharacters() = %v, want %v", gotChars, tt.wantChars)
			}
		})
	}
}

func Test_seqbag_ReverseComplement(t *testing.T) {
	type fields struct {
		seqmap          map[string]*seq
		seqs            []*seq
		ignoreidentical int
		alphabet        int
	}
	tests := []struct {
		name    string
		fields  fields
		wantseq []uint8
	}{
		{name: "t1",
			fields: fields{seqmap: nil,
				seqs:            []*seq{{sequence: []uint8("ATUGCYRSWKMBDHVN*.")}},
				ignoreidentical: IGNORE_NONE,
				alphabet:        NUCLEOTIDS},
			wantseq: []uint8(".*NBDHVKMWSYRGCAAT")},
		{name: "t2",
			fields: fields{seqmap: nil,
				seqs:            []*seq{{sequence: []uint8("TAACGRYSWMKVHDBN*.")}},
				ignoreidentical: IGNORE_NONE,
				alphabet:        NUCLEOTIDS},
			wantseq: []uint8(".*NVHDBMKWSRYCGTTA")},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			sb := &seqbag{
				seqmap:          tt.fields.seqmap,
				seqs:            tt.fields.seqs,
				ignoreidentical: tt.fields.ignoreidentical,
				alphabet:        tt.fields.alphabet,
			}
			sb.AutoAlphabet()
			sb.ReverseComplement()
			s, _ := sb.GetSequenceCharById(0)
			if !reflect.DeepEqual(s, tt.wantseq) {
				t.Errorf("seqbag.ReverseComplement() = %v, want %v", string(s), string(tt.wantseq))
			}
		})
	}
}
