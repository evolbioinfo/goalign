# Goalign: toolkit and api for alignment manipulation

## API

### append

Append alignments to an input alignment.

This commands adds the sequences of a set of alignments to a reference alignement
specified by -i.

If sequences do not have the same length than the reference alignment, then returns an error.

If format is phylip, it may contain several alignments in one file. 
Then we can append all of them at once:
goalign append -i refalign.phy aligns.phy

If format is Fasta, several alignments may be given in the form:
goalign append -i align.fasta others*.fasta

```go
package main

import (
	"bufio"
	"fmt"
	"io"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io/fasta"
	"github.com/evolbioinfo/goalign/io/utils"
)

func main() {
	var fi, fi2 io.Closer
	var r, r2 *bufio.Reader
	var err error
	var a1, a2 align.Alignment

	/* First alignment */
	/* Get reader (plain text or gzip) */
	if fi, r, err = utils.GetReader("align_1.fa"); err != nil {
		panic(err)
	}
	/* Parse Fasta */
	if a1, err = fasta.NewParser(r).Parse(); err != nil {
		panic(err)
	}
	fi.Close()

	/* Second alignment */
	/* Get reader (plain text or gzip) */
	if fi2, r2, err = utils.GetReader("align_2.fa"); err != nil {
		panic(err)
	}
	/* Parse Fasta */
	if a2, err = fasta.NewParser(r2).Parse(); err != nil {
		panic(err)
	}
	fi2.Close()

	/* Printing new alignment */
	if err = a1.Append(a2); err != nil {
		panic(err)
	}
	fmt.Println(fasta.WriteAlignment(a1))
}
```
