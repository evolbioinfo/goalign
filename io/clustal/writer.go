package clustal

import (
	"bytes"
	"fmt"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/version"
)

const (
	CLUSTAL_LINE = 50
)

func min_int(a int, b int) int {
	if a < b {
		return a
	}
	return b
}

func WriteAlignment(al align.Alignment) string {
	var buf bytes.Buffer
	cursize := 0

	// Get length of the longest name
	maxnamelength := 0
	al.IterateChar(func(name string, seq []rune) {
		if len(name) > maxnamelength {
			maxnamelength = len(name)
		}
	})

	buf.WriteString(fmt.Sprintf("CLUSTAL W (goalign version %s)\n\n", version.Version))
	for cursize < al.Length() {
		if cursize > 0 {
			buf.WriteRune('\n')
		}
		end := 0
		al.IterateChar(func(name string, seq []rune) {
			buf.WriteString(name)
			for i := len(name); i < maxnamelength+3; i++ {
				buf.WriteRune(' ')
			}

			end = min_int(cursize+CLUSTAL_LINE, len(seq))
			for j := cursize; j < end; j++ {
				buf.WriteRune(seq[j])
			}
			buf.WriteRune(' ')
			buf.WriteString(fmt.Sprintf("%d", end))
			buf.WriteRune('\n')
		})
		// Conservation line
		// White spaces
		for i := 0; i < maxnamelength+3; i++ {
			buf.WriteRune(' ')
		}
		// Each position in the line
		for pos := cursize; pos < end; pos++ {
			conservation, _ := al.SiteConservation(pos)
			switch conservation {
			case align.POSITION_IDENTICAL:
				buf.WriteRune('*')
			case align.POSITION_CONSERVED:
				buf.WriteRune(':')
			case align.POSITION_SEMI_CONSERVED:
				buf.WriteRune('.')
			default:
				buf.WriteRune(' ')
			}
		}

		buf.WriteRune('\n')
		cursize += CLUSTAL_LINE
	}
	return buf.String()
}
