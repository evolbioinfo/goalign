# Goalign: toolkit and api for alignment manipulation

## API

### unalign

Printing unaligned version of an input alignment

```go
package main

import (
	"bufio"
	"fmt"
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io"
	"github.com/fredericlemoine/goalign/io/fasta"
	"github.com/fredericlemoine/goalign/io/utils"
)

func main() {
	var fi *os.File
	var r *bufio.Reader
	var err error
	var al align.Alignment

	fi, r, err = utils.GetReader("align.fa")
	if err != nil {
		io.ExitWithMessage(err)
	}

	al, err = fasta.NewParser(r).Parse()
	if err != nil {
		io.ExitWithMessage(err)
	}
	fi.Close()

	fmt.Println(fasta.WriteSequences(al))
}
```
