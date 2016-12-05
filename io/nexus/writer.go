package nexus

import (
	"bytes"
	"fmt"
	"github.com/fredericlemoine/goalign/align"
)

func min_int(a int, b int) int {
	if a < b {
		return a
	}
	return b
}

func WriteAlignment(al align.Alignment) string {
	var buf bytes.Buffer

	var seqtype string = "nucleotide"

	if al.Alphabet() == align.AMINOACIDS {
		seqtype = "protein"
	}

	buf.WriteString("#NEXUS\n")
	buf.WriteString("begin data;\n")
	buf.WriteString(fmt.Sprintf("dimensions ntax=%d nchar=%d;\n", al.NbSequences(), al.Length()))
	buf.WriteString(fmt.Sprintf("format datatype=%s;\n", seqtype))
	buf.WriteString(fmt.Sprintf("matrix\n"))
	al.Iterate(func(name string, seq string) {
		buf.WriteString(name)
		buf.WriteString(" ")
		buf.WriteString(seq)
		buf.WriteRune('\n')
	})
	buf.WriteString(";\n")
	buf.WriteString("end;\n")

	return buf.String()
}
