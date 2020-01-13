package partition

type Token int64

var eof = rune(0)

const (
	ILLEGAL    Token = iota
	IDENTIFIER       // Identifier of model or partition
	SEPARATOR        // field separator : ,
	EQUAL            // Separator between model name and definition
	RANGE            // When defining a range ex 1-500
	MODULO           // Take one site every x sites: '/'
	DECIMAL          // Decimal
	ENDOFLINE        // End of line token
	EOF              // End of File
)

func isEndOfLine(ch rune) bool {
	return ch == '\n' || ch == '\r'
}

func isWhiteSpace(ch rune) bool {
	return ch == ' '
}

func isCR(ch rune) bool {
	return ch == '\r'
}

func isNL(ch rune) bool {
	return ch == '\n'
}

func isIdent(ch rune) bool {
	return ch != '\n' && ch != '\r' && ch != ',' && ch != '-' && ch != '/' && ch != '=' && ch != ' '
}
