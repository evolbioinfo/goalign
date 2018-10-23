# Goalign: toolkit and api for alignment manipulation

## API

### Orf

Detect the longest orf in an input dataset

```go
package main

import (
	"bufio"
	"fmt"
	"io"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io/fasta"
	"github.com/fredericlemoine/goalign/io/utils"
)

func main() {
	var fi io.Closer
	var r *bufio.Reader
	var err error
	var seqs align.SeqBag
	var orf align.Sequence

	// Get reader (plain text or gzip)
	if fi, r, err = utils.GetReader("seqs.fa"); err != nil {
		panic(err)
	}
	defer fi.Close()

	// Parse Fasta unaligned sequence file
	if seqs, err = fasta.NewParser(r).ParseUnalign(); err != nil {
		panic(err)
	}
	// Removing '-'
	seqs = seqs.Unalign()
	// Search for the longest orf
	if orf, err = seqs.LongestORF(); err != nil {
		panic(err)
	}
	// Print sequence
	fmt.Println(orf.Sequence())
}
```
