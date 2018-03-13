# Goalign: toolkit and api for alignment manipulation

## API

### trim name

Trimming names of sequences (to 10 characters)

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

	var namemap map[string]string

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

	/* Trim names */
	if namemap, err = al.TrimNames(10); err != nil {
		panic(err)
	} else {
		fmt.Println(fasta.WriteAlignment(al))
	}

	/* Old names => New names */
	for k, v := range namemap {
		fmt.Println(k + "=>" + v)
	}
}
```

Automatically generate names of sequences
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

	var namemap map[string]string

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

	/* Trim names */
	if namemap, err = al.TrimNamesAuto(); err != nil {
		panic(err)
	} else {
		fmt.Println(fasta.WriteAlignment(al))
	}

	/* Old names => New names */
	for k, v := range namemap {
		fmt.Println(k + "=>" + v)
	}
}
```

### trim seq

Trimming sequences (10 nucleotides from the right)

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

	/* Get reader (plain text or gzip) */
	fi, r, err = utils.GetReader("align.fa")
	if err != nil {
		panic(err)
	}

	/* Parse fasta */
	al, err = fasta.NewParser(r).Parse()
	if err != nil {
		panic(err)
	}
	fi.Close()

	/* Trim sequences */
	if err := al.TrimSequences(10, false); err != nil {
		panic(err)
	} else {
		fmt.Println(fasta.WriteAlignment(al))
	}
}
```
