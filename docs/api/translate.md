# Goalign: toolkit and api for alignment manipulation

## API

### translate

This command translates an input sequence into amino acids while removing first character.


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
	var al align.Alignment
	var translated align.Alignment

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

	/* Translate */
	if translated, err = al.Translate(1); err != nil {
		panic(err)
	} else {
		fmt.Println(fasta.WriteAlignment(translated))
	}
}
```
