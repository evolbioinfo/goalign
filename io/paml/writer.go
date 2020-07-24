package paml

import (
	"bytes"
	"fmt"

	"github.com/evolbioinfo/goalign/align"
)

const (
	PAML_LINE  = 60
	PAML_BLOCK = 10
)

func min_int(a int, b int) int {
	if a < b {
		return a
	}
	return b
}

func WriteAlignment(al align.Alignment) string {
	var buf bytes.Buffer

	buf.WriteString(fmt.Sprintf("  %d %d  I\n", al.NbSequences(), al.Length()))
	al.Iterate(func(name string, seq string) bool {
		buf.WriteString(name + "\n")
		return false
	})
	buf.WriteRune('\n')
	cursize := 0
	for cursize < al.Length() {
		if cursize > 0 {
			buf.WriteString(fmt.Sprintf("%d\n", cursize+1))
		}
		al.IterateChar(func(name string, seq []uint8) bool {
			for i := cursize; i < cursize+PAML_LINE && i < len(seq); i += PAML_BLOCK {
				if i > cursize {
					buf.WriteString(" ")
				}
				end := min_int(i+PAML_BLOCK, len(seq))
				for j := i; j < end; j++ {
					buf.WriteByte(seq[j])
				}
			}
			buf.WriteString("\n")
			return false
		})
		cursize += PAML_LINE
	}
	return buf.String()
}
