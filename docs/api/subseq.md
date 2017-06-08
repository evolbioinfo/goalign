# Goalign: toolkit and api for alignment manipulation

## API

### subseq

Extracting sub-alignment (position 10 0-based inclusive, and with length 15) from an input alignment:

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
	var subalign align.Alignment = nil

	/* Get reader (plain text or gzip) */
	fi, r, err = utils.GetReader("align.fa")
	if err != nil {
		io.ExitWithMessage(err)
	}

	/* Parse Fasta */
	al, err = fasta.NewParser(r).Parse()
	if err != nil {
		io.ExitWithMessage(err)
	}
	fi.Close()

	/* Subalignment from position 10 (0-based inclusive),
	and with length 15 */
	subalign, err = al.SubAlign(10, 15)
	if err != nil {
		io.ExitWithMessage(err)
	}

	fmt.Print(fasta.WriteAlignment(subalign))
}
```
