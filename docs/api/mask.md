# Goalign: toolkit and api for alignment manipulation

## API

### mask

Mask part of an input fasta alignment

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

	/* First Alignment*/

	/* Get reader (plain text or gzip) */
	if fi, r, err = utils.GetReader("align.fa"); err != nil {
		panic(err)
	}
	defer fi.Close()

	/* Parse Fasta */
	if al, err = fasta.NewParser(r).Parse(); err != nil {
		panic(err)
	}

	if err = al.Mask(0, 2); err != nil {
		panic(err)
	}

	fmt.Println(fasta.WriteAlignment(al))
}
```
