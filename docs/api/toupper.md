# Goalign: toolkit and api for alignment manipulation

## API

### To Upper case

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
	var alignmentReader *bufio.Reader
	var err error
	var al align.Alignment

	/* Parse Fasta */
	if fi, alignmentReader, err = utils.GetReader("align.fa"); err != nil {
		panic(err)
	}
	defer fi.Close()

	if al, err = fasta.NewParser(alignmentReader).Parse(); err != nil {
		panic(err)
	}
	defer fi.Close()

    al.ToUpper()

    fmt.Println(fasta.WriteAlignment(al))
}
```
