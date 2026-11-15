# Goalign: toolkit and api for alignment manipulation

## API

### compute entropy

`Entropy` and `EntropyFiltered` return the Shannon entropy of one site, in **nats** (natural logarithm). Divide by `math.Log(2)` to get bits.

- `Entropy(site, removegaps)`: `*` and `.` are never counted, and `-` is not counted if `removegaps` is true. It is the same as `EntropyFiltered(site, removegaps, false)`.
- `EntropyFiltered(site, removegaps, removeambiguities)`: if `removeambiguities` is true, it only counts standard characters, depending on the alignment alphabet: `A C G T U` for nucleotides, and the 20 standard amino acids for proteins. N, X, IUPAC ambiguity codes, B, Z, J, etc. are not counted. It returns an error if the alignment alphabet is unknown.

Characters are upper-cased before counting. A site where no character is counted has an entropy of `NaN`. A site outside the alignment (`< 0` or `>= Length()`) returns an error.

```go
package main

import (
	"bufio"
	"fmt"
	"io"
	"math"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io/fasta"
	"github.com/evolbioinfo/goalign/io/utils"
)

func main() {
	var fi io.Closer
	var r *bufio.Reader
	var err error
	var al align.Alignment
	var e float64

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

	/* Entropy of each site, in bits, ignoring gaps and ambiguities */
	for site := 0; site < al.Length(); site++ {
		if e, err = al.EntropyFiltered(site, true, true); err != nil {
			panic(err)
		}
		fmt.Printf("%d\t%.6f\n", site, e/math.Log(2))
	}
}
```

### compute pssm

`Pssm(log, pseudocount, normalization)` computes a position specific scoring matrix, as a map from character to per-site values. Possible values of `normalization` are `align.PSSM_NORM_NONE`, `align.PSSM_NORM_FREQ`, `align.PSSM_NORM_DATA`, `align.PSSM_NORM_UNIF` and `align.PSSM_NORM_LOGO`.

With `align.PSSM_NORM_LOGO`, the site entropy is computed in **bits** (log2), as is usual for sequence logos. `Entropy` and `EntropyFiltered`, in contrast, return nats.

```go
	var pssm map[uint8][]float64
	if pssm, err = al.Pssm(false, 0.0, align.PSSM_NORM_LOGO); err != nil {
		panic(err)
	}
	for c, values := range pssm {
		fmt.Printf("%c\t%v\n", c, values)
	}
```
