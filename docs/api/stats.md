# Goalign: toolkit and api for alignment manipulation

## API

### stats

Printing statistics about an input alignment:

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

	/* Print alignment length*/
	fmt.Printf("Length=%d\n", al.Length())
	/* Print number of sequences */
	fmt.Printf("#Seqs=%d\n", al.NbSequences())
	/* Print avg allements per sites */
	fmt.Printf("Avg alleles/sites=%f\n", al.AvgAllelesPerSite())
	/* Print number of occurences of each characters */
	for nt, nb := range al.CharStats() {
		fmt.Printf("%c = %d\n", nt, nb)
	}
}
```
