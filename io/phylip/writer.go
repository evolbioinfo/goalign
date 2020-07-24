package phylip

import (
	"bytes"
	"fmt"

	"github.com/evolbioinfo/goalign/align"
)

const (
	PHYLIP_LINE  = 60
	PHYLIP_BLOCK = 10
)

func min_int(a int, b int) int {
	if a < b {
		return a
	}
	return b
}

func WriteAlignment(al align.Alignment, strict, oneline, noblock bool) string {
	var buf bytes.Buffer
	var header bool = true
	var line_length = PHYLIP_LINE
	var block_length = PHYLIP_BLOCK

	cursize := 0
	buf.WriteString(fmt.Sprintf("   %d   %d\n", al.NbSequences(), al.Length()))

	if oneline {
		line_length = al.Length()
	}
	if noblock {
		block_length = line_length
	}

	for cursize < al.Length() {
		if cursize > 0 {
			buf.WriteString("\n")
		}
		al.IterateChar(func(name string, seq []uint8) bool {
			if header {
				if strict {
					buf.WriteString(fmt.Sprintf("%-10s", name[:min_int(10, len(name))]))
				} else {
					buf.WriteString(name)
					buf.WriteString("  ")
				}
			}

			for i := cursize; i < cursize+line_length && i < len(seq); i += block_length {
				if i > cursize {
					buf.WriteString(" ")
				} else if !header {
					if strict {
						buf.WriteString("          ")
					} else {
						buf.WriteString("   ")
					}
				}
				end := min_int(i+block_length, len(seq))
				for j := i; j < end; j++ {
					buf.WriteByte(seq[j])
				}
			}
			buf.WriteString("\n")
			return false
		})
		cursize += line_length
		header = false
	}
	return buf.String()
}
