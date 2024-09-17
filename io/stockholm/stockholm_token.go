package stockholm

type Token int64

var eof = rune(0)

const (
	ILLEGAL Token = iota
	EOF
	WS
	IDENT     // Name of Node, or comment, or keyword
	NUMERIC   // Any numerical value
	ENDOFLINE // \r \n

	// Keywords
	STOCKHOLM // STOCKHOLM  : Start of Stockholm file
	MARKUP    // #
	END       // //
	TREE      // A specific tree in the BEGIN TREES section
)

func isWhitespace(ch rune) bool {
	return ch == ' ' || ch == '\t'
}

func isIdent(ch rune) bool {
	return ch != '[' && ch != ']' && ch != ';' && ch != '=' && ch != '\r' && ch != '\n' && !isWhitespace(ch)
}

func isEndOfLine(ch rune) bool {
	return ch == '\n' || ch == '\r'
}

func isCR(ch rune) bool {
	return ch == '\r'
}

func isNL(ch rune) bool {
	return ch == '\n'
}
