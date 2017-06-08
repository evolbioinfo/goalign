# Goalign: toolkit and api for alignment manipulation

## API

### trim name

Trimming names of sequences

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

	var namemap map[string]string

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

	/* Trim names */
	if namemap, err = al.TrimNames(3); err != nil {
		io.ExitWithMessage(err)
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

	/* Parse fasta */
	al, err = fasta.NewParser(r).Parse()
	if err != nil {
		io.ExitWithMessage(err)
	}
	fi.Close()

	/* Trim sequences */
	if err := al.TrimSequences(10, false); err != nil {
		io.ExitWithMessage(err)
	} else {
		fmt.Println(fasta.WriteAlignment(al))
	}
}
```
