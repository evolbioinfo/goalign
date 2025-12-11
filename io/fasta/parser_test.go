package fasta

import (
	"bytes"
	"fmt"
	"strings"
	"testing"
)

var fastastring string = ">s1\nACGATCGATTACTACTGAC\nACGACTGATCGATCG"
var fastastring2 string = " >s1\nACGATCGATTACTACTGAC\nACGACTGATCGATCG"
var fastastring3 string = ">s1\nACGATCGATTACTACTGAC\nACGACTGATCGATCG\n>s2\nACGATCGATTACTACTGAC\nACGACTGATCGATCG\n"
var fastastring4 string = ">s1\nACGATCGATTACTACTGAC\nACGACTGATCGATCG\n>s2\nACGATCGATTACTACTGAC\nACGACTGATCGATC\n"
var seq = "AACGTACGTACAGCTAGCTATGTACTGATCATGCTAGCTGC\nACCAGCATGCTACTACTAGCTCGATGCATCGCATATGCAC\n"

func TestParse(t *testing.T) {
	align, err := NewParser(strings.NewReader(fastastring)).Parse()

	if err != nil {
		t.Error(err)
	}
	if align.Length() != 34 {
		t.Errorf("Alignment length is not 34 %d", align.Length())
	}
	if align.NbSequences() != 1 {
		t.Errorf("There is not 1 sequence in the alignment %d", align.NbSequences())
	}

	_, err2 := NewParser(strings.NewReader(fastastring2)).Parse()

	if err2 == nil {
		t.Errorf("There should be an error while parsing fastastring2")
	}

	align3, err3 := NewParser(strings.NewReader(fastastring3)).Parse()

	if err3 != nil {
		t.Error(err3)
	}

	if align3.NbSequences() != 2 {
		t.Errorf("There are not 2 sequence in the alignment %d", align3.NbSequences())
	}

	_, err4 := NewParser(strings.NewReader(fastastring4)).Parse()

	if err4 == nil {
		t.Errorf("There should be an error while parsing fastastring4, which has different length sequences")
	}

	var fasta bytes.Buffer
	for i := 0; i < 1000; i++ {
		fasta.WriteString(fmt.Sprintf(">s%d\n%s", i, seq))
	}
	align5, err5 := NewParser(strings.NewReader(fasta.String())).Parse()
	if err5 != nil {
		t.Error(err5)
	}
	if align5.NbSequences() != 1000 {
		t.Errorf("Alignment has not 1000 sequences : %d", align5.NbSequences())
	}
}
