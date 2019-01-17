# Goalign: toolkit and api for alignment manipulation

## API

### replace

Replace characters in sequences of input alignment using a regex.

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
	var al align.SeqBag
	// or var al align.Alignment if aligned sequences

	/* First Alignment*/

	/* Get reader (plain text or gzip) */
	if fi, r, err = utils.GetReader("align.fa"); err != nil {
		panic(err)
	}
	defer fi.Close()

	/* Parse Fasta */
	if al, err = fasta.NewParser(r).ParseUnalign(); err != nil {
	/* for aligned seqs: if al, err = fasta.NewParser(r).Parse(); err != nil {*/
		panic(err)
	}

	if err = al.Replace("GA.", "---", true); err != nil {
		panic(err)
	}

	fmt.Println(fasta.WriteAlignment(al))
}
```
