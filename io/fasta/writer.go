package fasta

import (
	"bytes"
	"github.com/fredericlemoine/goalign/align"
)

const (
	FASTA_LINE = 80
)

func min_int(a int, b int) int {
	if a < b {
		return a
	}
	return b
}

func WriteAlignment(al align.Alignment) string {
	var buf bytes.Buffer

	al.IterateChar(func(name string, seq []rune) {
		buf.WriteString(">")
		buf.WriteString(name)
		buf.WriteString("\n")
		for i := 0; i < len(seq); i += FASTA_LINE {
			end := min_int(i+FASTA_LINE, len(seq))
			for j := i; j < end; j++ {
				buf.WriteRune(seq[j])
			}
			buf.WriteString("\n")
		}
	})
	return buf.String()
}
