# Goalign: toolkit and api for alignment manipulation

## API

### reformat

Reformat a fasta input alignment into different formats

```go
package main

import (
	"bufio"
	"fmt"
	"io"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io/fasta"
	"github.com/fredericlemoine/goalign/io/nexus"
	"github.com/fredericlemoine/goalign/io/paml"
	"github.com/fredericlemoine/goalign/io/phylip"
	"github.com/fredericlemoine/goalign/io/clustal"
	"github.com/fredericlemoine/goalign/io/utils"
)

func main() {
	var fi io.Closer
	var r *bufio.Reader
	var err error
	var al align.Alignment

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

	/* Printing FASTA alignment */
	fmt.Println(fasta.WriteAlignment(al))
	/* Printing PHYLIP alignment */
	fmt.Println(phylip.WriteAlignment(al, false))
	/* Printing NEXUS alignment */
	fmt.Println(nexus.WriteAlignment(al))
	/* Printing PAML format alignment */
	fmt.Println(paml.WriteAlignment(al))
    /* Printing Clustal format alignment */
	fmt.Println(clustal.WriteAlignment(al))
}
```
