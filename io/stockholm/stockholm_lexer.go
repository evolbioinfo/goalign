package stockholm

import (
	"bufio"
	"bytes"
	"io"
	"strconv"
	"strings"

	aio "github.com/evolbioinfo/goalign/io"
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

	// If we see whitespace then consume all contiguous whitespace.
	// If we see a letter then consume as an ident or reserved word.
	if isWhitespace(ch) {
		s.unread()
		return s.scanWhitespace()
	}

	if isEndOfLine(ch) {
		if isCR(ch) {
			ch := s.read()
			if isNL(ch) {
				return ENDOFLINE, ""
			} else {
				aio.PrintMessage("\\r without \\n detected...")
			}
		} else {
			return ENDOFLINE, ""
		}
	}

	switch ch {
	case eof:
		return EOF, ""
	case '#':
		return MARKUP, string(ch)
	}

	s.unread()
	return s.scanIdent()
}

// scanWhitespace consumes the current rune and all contiguous whitespace.
func (s *Scanner) scanWhitespace() (tok Token, lit string) {
	// Create a buffer and read the current character into it.
	var buf bytes.Buffer
	buf.WriteRune(s.read())

	// Read every subsequent whitespace character into the buffer.
	// Non-whitespace characters and EOF will cause the loop to exit.
	for {
		if ch := s.read(); ch == eof {
			break
		} else if !isWhitespace(ch) {
			s.unread()
			break
		} else {
			buf.WriteRune(ch)
		}
	}

	return WS, buf.String()
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

	_, err := strconv.ParseInt(buf.String(), 10, 64)
	if err != nil {
		switch strings.ToUpper(buf.String()) {
		case "STOCKHOLM":
			return STOCKHOLM, buf.String()
		case "//":
			return END, buf.String()
		default:
			return IDENT, buf.String()
		}
	} else {
		return NUMERIC, buf.String()
	}
}
