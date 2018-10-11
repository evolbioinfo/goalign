# Goalign: toolkit and api for alignment manipulation

## API

### identical

Telling whether two alignments are identical

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
	var fi, fi2 io.Closer
	var r, r2 *bufio.Reader
	var err error
	var al, al2 align.Alignment

	/* First Alignment*/

	/* Get reader (plain text or gzip) */
	if fi, r, err = utils.GetReader("align.fa"); err != nil {
		panic(err)
	}

	/* Parse Fasta */
	if al, err = fasta.NewParser(r).Parse(); err != nil {
		panic(err)
	}
	fi.Close()

	/* Second Alignments */

	/* Get reader (plain text or gzip) */
	if fi2, r2, err = utils.GetReader("align2.fa"); err != nil {
		panic(err)
	}

	/* Parse Fasta */
	if al2, err = fasta.NewParser(r2).Parse(); err != nil {
		panic(err)
	}
	fi2.Close()

	/* Printing unaligned sequences */
	fmt.Println(al.Identical(al2))
}
```
