package stockholm

import (
	"bytes"

	"github.com/evolbioinfo/goalign/align"
)

func WriteAlignment(al align.Alignment) string {
	var buf bytes.Buffer

	buf.WriteString("# STOCKHOLM 1.0\n")
	buf.WriteString("#=GF ID   Goalign generated alignment\n")
	al.Iterate(func(name string, seq string) bool {
		buf.WriteString(name)
		buf.WriteString("\t")
		buf.WriteString(seq)
		buf.WriteRune('\n')
		return false
	})
	buf.WriteString("//")

	return buf.String()
}
