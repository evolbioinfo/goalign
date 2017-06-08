# Goalign: toolkit and api for alignment manipulation

## API

### subset

Extracting sequences ""Seq0001" and "Seq0002" from an input alignment

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
	var filtered align.Alignment = nil
	var subset map[string]bool

	/* Sequence names to keep */
	subset = make(map[string]bool)
	for _, name := range []string{"Seq0001", "Seq0002"} {
		subset[name] = true
	}

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

	/* Iterate over alignment sequences */
	al.Iterate(func(name string, sequence string) {
		if filtered == nil {
			filtered = align.NewAlign(align.DetectAlphabet(sequence))
		}
		/* Adding only the desired one to the filtered alignment */
		if _, ok := subset[name]; ok {
			filtered.AddSequence(name, sequence, "")
		}
	})

	fmt.Print(fasta.WriteAlignment(filtered))
}
```
