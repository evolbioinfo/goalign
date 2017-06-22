package phylip

import (
	"fmt"
	"strings"
	"testing"
)

var phylipstring1 string = "2 20\nseq1 ACGATCGATC\nseq2 AGCTCGTCGA\n\nACGCTAGCTA\nGATCGCTAGC\n"
var phylipstring2 string = "2 20\nseq1      ACGATCGATC\nseq1234567AGCTCGTCGA\n\nACGCTAGCTA\nGATCGCTAGC\n"

func TestParse(t *testing.T) {
	align, err := NewParser(strings.NewReader(phylipstring1), false).Parse()

	if err != nil {
		t.Error(err)
	}
	if align.Length() != 20 {
		t.Error("Alignment length is not 20" + fmt.Sprintf("%d", align.Length()))
	}
	if align.NbSequences() != 2 {
		t.Error("There is not 2 sequence in the alignment" + fmt.Sprintf("%d", align.NbSequences()))
	}

	// _, err2 := NewParser(strings.NewReader(fastastring2)).Parse()

	// if err2 == nil {
	// 	t.Error("There should be an error while parsing fastastring2")
	// }

	// align3, err3 := NewParser(strings.NewReader(fastastring3)).Parse()

	// if err3 != nil {
	// 	t.Error(err3)
	// }

	// if align3.NbSequences() != 2 {
	// 	t.Error("There are not 2 sequence in the alignment" + fmt.Sprintf("%d", align3.NbSequences()))
	// }

	// _, err4 := NewParser(strings.NewReader(fastastring4)).Parse()

	// if err4 == nil {
	// 	t.Error("There should be an error while parsing fastastring4, which has different length sequences")
	// }

	// var fasta bytes.Buffer
	// for i := 0; i < 1000; i++ {
	// 	fasta.WriteString(">s" + fmt.Sprint("%d", i) + "\n")
	// 	fasta.WriteString(seq)
	// }
	// align5, err5 := NewParser(strings.NewReader(fasta.String())).Parse()
	// if err5 != nil {
	// 	t.Error(err5)
	// }
	// if align5.NbSequences() != 1000 {
	// 	t.Error("Alignment has not 1000 sequences : " + fmt.Sprintf("%d", align5.NbSequences()))
	// }
}

func TestParseStrict(t *testing.T) {
	align, err := NewParser(strings.NewReader(phylipstring2), true).Parse()

	if err != nil {
		t.Error(err)
	}
	if align.Length() != 20 {
		t.Error("Alignment length is not 20" + fmt.Sprintf("%d", align.Length()))
	}
	if align.NbSequences() != 2 {
		t.Error("There is not 2 sequence in the alignment" + fmt.Sprintf("%d", align.NbSequences()))
	}

	// _, err2 := NewParser(strings.NewReader(fastastring2)).Parse()

	// if err2 == nil {
	// 	t.Error("There should be an error while parsing fastastring2")
	// }

	// align3, err3 := NewParser(strings.NewReader(fastastring3)).Parse()

	// if err3 != nil {
	// 	t.Error(err3)
	// }

	// if align3.NbSequences() != 2 {
	// 	t.Error("There are not 2 sequence in the alignment" + fmt.Sprintf("%d", align3.NbSequences()))
	// }

	// _, err4 := NewParser(strings.NewReader(fastastring4)).Parse()

	// if err4 == nil {
	// 	t.Error("There should be an error while parsing fastastring4, which has different length sequences")
	// }

	// var fasta bytes.Buffer
	// for i := 0; i < 1000; i++ {
	// 	fasta.WriteString(">s" + fmt.Sprint("%d", i) + "\n")
	// 	fasta.WriteString(seq)
	// }
	// align5, err5 := NewParser(strings.NewReader(fasta.String())).Parse()
	// if err5 != nil {
	// 	t.Error(err5)
	// }
	// if align5.NbSequences() != 1000 {
	// 	t.Error("Alignment has not 1000 sequences : " + fmt.Sprintf("%d", align5.NbSequences()))
	// }
}
