# Goalign: toolkit and api for alignment manipulation

## API

### subseq

Extracting sub-alignment (position 10 0-based inclusive, and with length 15) from an input alignment:

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
	var fi io.Closer
	var r *bufio.Reader
	var err error
	var al align.Alignment
	var subalign align.Alignment = nil

	/* Get reader (plain text or gzip) */
	fi, r, err = utils.GetReader("align.fa")
	if err != nil {
		panic(err)
	}

	/* Parse Fasta */
	al, err = fasta.NewParser(r).Parse()
	if err != nil {
		panic(err)
	}
	fi.Close()

	/* Subalignment from position 10 (0-based inclusive),
	and with length 15 */
	subalign, err = al.SubAlign(10, 15)
	if err != nil {
		panic(err)
	}

	/* Printing alignment in Fasta */
	fmt.Print(fasta.WriteAlignment(subalign))
}
```
