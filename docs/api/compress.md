# Goalign: toolkit and api for alignment manipulation

## API

### compress

Remove identical patterns/sites

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
	var weights []int

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

	/* Compress */
	weights = al.Compress()
	fmt.Println(fasta.WriteAlignment(al))
	for _, w := range weights {
		fmt.Println(w)
	}
}
```
