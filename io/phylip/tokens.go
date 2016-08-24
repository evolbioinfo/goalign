package phylip

type Token int64

var eof = rune(0)

const (
	ILLEGAL    Token = iota
	INTEGER          // Number of sequences or length
	IDENTIFIER       // Identifier of sequence or part of sequence
	ENDOFLINE        // End of line token
	EOF              // End of File
	WS               // Whitespace
	NUMERIC          // Number of taxa and length of sequences
)

func isEndOfLine(ch rune) bool {
	return ch == '\n' || ch == '\r'
}

func isIdent(ch rune) bool {
	return ch != '\n' && ch != ' ' && ch != '\n' && ch != '\r'
}

func isWhitespace(ch rune) bool {
	return ch == ' ' || ch == '\t'
}
