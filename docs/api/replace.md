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

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io/fasta"
	"github.com/evolbioinfo/goalign/io/utils"
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

Replace characters in sequences of input alignment using its position+sequence name.

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
	var al align.SeqBag
	// or var al align.Alignment if aligned sequences

	/* First Alignment*/

	/* Get reader (plain text or gzip) */
	if fi, r, err = utils.GetReader("align.fa"); err != nil {
		panic(err)
	}
	defer fi.Close()

	/* Parse Fasta */
	if al, err = fasta.NewParser(r).Parse(); err != nil {
	/* for aligned seqs: if al, err = fasta.NewParser(r).Parse(); err != nil {*/
		panic(err)
	}

	// Will write a G at position 9 (0 based) of sequence "Seq001"
	if err = al.ReplaceChar("Seq001", 9, "G"); err != nil {
		panic(err)
	}

	fmt.Println(fasta.WriteAlignment(al))
}
```
