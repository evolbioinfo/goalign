# Goalign: toolkit and api for alignment manipulation

## API

### splits

Printing statistics about an input alignment:

```go
package main

import (
	"bufio"
	"fmt"
	"io"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io/fasta"
	"github.com/evolbioinfo/goalign/io/partition"
	"github.com/evolbioinfo/goalign/io/utils"
)

func main() {
	var fi io.Closer
	var partitionReader *bufio.Reader
	var alignmentReader *bufio.Reader
	var partitionParser *partition.Parser
	var splitAligns []align.Alignment
	var ps *align.PartitionSet
	var err error
	var al align.Alignment

	/* Parse Fasta */
	if fi, alignmentReader, err = utils.GetReader("align.fa"); err != nil {
		panic(err)
	}
	defer fi.Close()

	if al, err = fasta.NewParser(alignmentReader).Parse(); err != nil {
		panic(err)
	}
	defer fi.Close()

	/* Parse Partitions */
	if fi, partitionReader, err = utils.GetReader("partitions.txt"); err != nil {
		panic(err)
	}
	defer fi.Close()
	partitionParser = partition.NewParser(partitionReader)
	if ps, err = partitionParser.Parse(al.Length()); err != nil {
		panic(err)
	}
	if err = ps.CheckSites(); err != nil {
		panic(err)
	}

	/*  Split alignment per partition and write them on std out */
	if splitAligns, err = al.Split(ps); err != nil {
		panic(err)
	}

	for i, a := range splitAligns {
		fmt.Printf("Alignment %d\n", i)
		fmt.Println(fasta.WriteAlignment(a))
	}
}
```
