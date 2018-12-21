# Goalign: toolkit and api for alignment manipulation

## API

### phase

Printing amino acid sequence "phased" version of an input file

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
	var fi, fi2 io.Closer
	var r, r2 *bufio.Reader
	var err error
	var seqs, refs align.SeqBag
	var phaser align.Phaser
	var phased chan align.PhasedSequence

	/* Get ref sequences reader (plain text or gzip) */
	if fi, r, err = utils.GetReader("refs.fa"); err != nil {
		panic(err)
	}
	defer fi.Close()

	/* Parse Reference orfs (nt or aa) */
	if refs, err = fasta.NewParser(r).ParseUnalign(); err != nil {
		panic(err)
	}

	/* Get sequences reader (plain text or gzip) */
	if fi2, r2, err = utils.GetReader("align.fa"); err != nil {
		panic(err)
	}
	defer fi2.Close()

	/* Parse Fasta sequences */
	seqs, err = fasta.NewParser(r2).ParseUnalign()
	if err != nil {
		panic(err)
	}

	/* Initialize phaser */
	phaser = align.NewPhaser()
	phaser.SetLenCutoff(0.8)
	phaser.SetMatchCutoff(0.8)
	phaser.SetReverse(true)
	phaser.SetCutEnd(true)
	phaser.SetCpus(1)
	phaser.SetTranslate(true)

	/* Phase */
	if phased, err = phaser.Phase(refs, seqs); err != nil {
		panic(err)
	}

	for p := range phased {
		if p.Err != nil {
			panic(err)
		}
		if p.Removed {
			fmt.Printf("%s: Removed\n", p.NtSeq.Name())
		} else {
			fmt.Printf("Nt Sequence %s : %s\n", p.NtSeq.Name(), p.NtSeq.Sequence())
			fmt.Printf("Aa Sequence %s : %s\n", p.AaSeq.Name(), p.AaSeq.Sequence())
			fmt.Printf("Position: %d\n", p.Position)
			fmt.Printf("Length: %d\n", p.AaSeq.Length())
		}
	}
}
```
