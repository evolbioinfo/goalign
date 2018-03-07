# Goalign: toolkit and api for alignment manipulation

## API

### draw

Drawing an alignment in an HTML file using [BioJS](http://msa.biojs.net/).

```go
package main

import (
	"bufio"
	"io"
	"os"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/draw"
	"github.com/fredericlemoine/goalign/io/fasta"
	"github.com/fredericlemoine/goalign/io/utils"
)

func main() {
	var infile io.Closer
	var outfile *os.File

	var r *bufio.Reader
	var err error

	var al align.Alignment
	var l draw.AlignLayout

	/* Get reader (plain text or gzip) */
	infile, r, err = utils.GetReader("align.fa")
	if err != nil {
		panic(err)
	}

	/* Parse Fasta */
	al, err = fasta.NewParser(r).Parse()
	if err != nil {
		panic(err)
	}
	infile.Close()

	outfile, err = os.Create("align.html")
	w := bufio.NewWriter(outfile)
	l = draw.NewBioJSLayout(w)
	l.DrawAlign(al)
	w.Flush()
	outfile.Close()
}
```
