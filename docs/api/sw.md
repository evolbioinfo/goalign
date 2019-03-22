# Goalign: toolkit and api for alignment manipulation

## API

### Smith & Waterman

```go
package main

import (
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io/fasta"
)

func main() {
	seq1 := align.NewSequence("seq1", []rune("CTGGGGTTTAACCAGCCATGCCAGTGCAGGTTTAAGAACCGATCCGTACTCTGGGTTACTGATGAAGGATGGGCCGTATCGCCCCCTTGCGACGTTTCCA"), "")
	seq2 := align.NewSequence("seq2", []rune("TATTATCGTATCGTTTGCATAGACCCGTTATGCCAGCAGATACAGCGTCACAAACTTAGGCTGTAGGGCGTTAGCGGCGCTCCATGTTTAGACTCACGCC"), "")
	aligner := align.NewPwAligner(seq1, seq2, align.ALIGN_ALGO_SW)
	aligner.SetGapOpenScore(-10.0)
	aligner.SetGapExtendScore(-0.5)
	if al, err := aligner.Alignment(); err != nil {
		panic(err)
	} else {
		fmt.Println(fasta.WriteAlignment(al))
	}
}
```
