# Goalign: toolkit and api for alignment manipulation

## API

### phase

Printing ATG "phased" version of an input file

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
	var seqs, phased, aaphased align.SeqBag

	/* Get reader (plain text or gzip) */
	fi, r, err = utils.GetReader("align.fa")
	if err != nil {
		panic(err)
	}

	/* Parse Fasta */
	seqs, err = fasta.NewParser(r).ParseUnalign()
	if err != nil {
		panic(err)
	}
	fi.Close()

	if phased, aaphased, _, err = seqs.Phase(nil); err != nil {
		panic(err)
	}

	/* Printing nt unaligned phased sequences */
	fmt.Println(fasta.WriteAlignment(phased))

	/* Printing aa unaligned phased sequences */
	fmt.Println(fasta.WriteAlignment(aaphased))
}
```
