package fasta

import (
	"bytes"

	"github.com/evolbioinfo/goalign/align"
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

func WriteAlignment(sb align.SeqBag) string {
	var buf bytes.Buffer
	sb.IterateChar(func(name string, seq []rune) bool {
		buf.WriteString(">")
		buf.WriteString(name)
		buf.WriteString("\n")
		for i := 0; i < len(seq); i++ {
			if i%FASTA_LINE == 0 && i > 0 {
				buf.WriteString("\n")
			}
			buf.WriteRune(seq[i])
		}
		buf.WriteRune('\n')
		return false
	})
	return buf.String()
}

// Write input alignment as standard fasta sequences
// It removes "-" characters.
func WriteSequences(sb align.SeqBag) string {
	var buf bytes.Buffer

	sb.IterateChar(func(name string, seq []rune) bool {
		buf.WriteString(">")
		buf.WriteString(name)
		buf.WriteString("\n")
		nbchar := 0
		for i := 0; i < len(seq); i++ {
			if seq[i] != '-' {
				buf.WriteRune(seq[i])
				nbchar++
				if nbchar == FASTA_LINE {
					buf.WriteString("\n")
					nbchar = 0
				}
			}
		}
		if nbchar != 0 {
			buf.WriteString("\n")
		}
		return false
	})
	return buf.String()
}
