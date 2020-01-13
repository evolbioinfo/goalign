package partition

import (
	"bufio"
	"bytes"
	"io"
	"strconv"
)

// Scanner represents a lexical scanner.
type Scanner struct {
	r *bufio.Reader
}

// NewScanner returns a new instance of Scanner.
func NewScanner(r io.Reader) *Scanner {
	return &Scanner{r: bufio.NewReader(r)}
}

// read reads the next rune from the bufferred reader.
// Returns the rune(0) if an error occurs (or io.EOF is returned).
func (s *Scanner) read() rune {
	ch, _, err := s.r.ReadRune()
	if err != nil {
		return eof
	}
	return ch
}

// unread places the previously read rune back on the reader.
func (s *Scanner) unread() {
	_ = s.r.UnreadRune()
}

// Scan returns the next token and literal value.
func (s *Scanner) Scan() (tok Token, lit string) {
	// Read the next rune.
	ch := s.read()

	if isEndOfLine(ch) {
		if isCR(ch) {
			ch := s.read()
			if !isNL(ch) {
				s.unread()
			}
		}
		return ENDOFLINE, ""
	}

	for isWhiteSpace(ch) {
		ch = s.read()
	}

	switch ch {
	case eof:
		return EOF, ""
	case ',':
		return SEPARATOR, string(ch)
	case '=':
		return EQUAL, string(ch)
	case '-':
		return RANGE, string(ch)
	case '/':
		return MODULO, string(ch)
	}

	s.unread()
	tok, ident := s.scanIdent()

	_, err := strconv.ParseInt(ident, 10, 64)
	if err != nil {
		return IDENTIFIER, ident
	} else {
		return DECIMAL, ident
	}
}

// scanIdent consumes the current rune and all contiguous ident runes.
func (s *Scanner) scanIdent() (tok Token, lit string) {
	// Create a buffer and read the current character into it.
	var buf bytes.Buffer
	buf.WriteRune(s.read())

	// Read every subsequent ident character into the buffer.
	// Non-ident characters and EOF will cause the loop to exit.
	for {
		if ch := s.read(); ch == eof {
			break
		} else if !isIdent(ch) {
			s.unread()
			break
		} else {
			_, _ = buf.WriteRune(ch)
		}
	}

	return IDENTIFIER, buf.String()
}
