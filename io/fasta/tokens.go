package fasta

type Token int64

var eof = rune(0)

const (
	ILLEGAL    Token = iota
	STARTIDENT       // > start of ident line
	IDENTIFIER       // Identifier of sequence or part of sequence
	ENDOFLINE        // End of line token
	EOF              // End of File
)

func isEndOfLine(ch rune) bool {
	return ch == '\n' || ch == '\r'
}

func isIdent(ch rune) bool {
	return ch != '\n' && ch != '\r'
}
