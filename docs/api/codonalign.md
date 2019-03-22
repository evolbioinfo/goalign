# Goalign: toolkit and api for alignment manipulation

## API

### codonalign

Aligns a given nt fasta file using a corresponding aa alignment.

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
	var aa, codonaligned align.Alignment
	var nt align.SeqBag

	/* Amino Acid Alignment */

	/* Get reader (plain text or gzip) */
	if fi, r, err = utils.GetReader("align_aa.fa"); err != nil {
		panic(err)
	}

	/* Parse Fasta */
	if aa, err = fasta.NewParser(r).Parse(); err != nil {
		panic(err)
	}
	fi.Close()

	/* Nt fasta sequences */

	/* Get reader (plain text or gzip) */
	if fi2, r2, err = utils.GetReader("align_nt.fa"); err != nil {
		panic(err)
	}

	/* Parse Fasta */
	if nt, err = fasta.NewParser(r2).ParseUnalign(); err != nil {
		panic(err)
	}
	fi2.Close()

	/* Printing codon alingned nt sequences */
	if codonaligned, err = aa.CodonAlign(nt); err != nil {
		panic(err)
	}
	fmt.Println(fasta.WriteAlignment(codonaligned))
}
```
