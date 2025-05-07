package clustal

import (
	"bufio"
	"bytes"
	"errors"
	"io"
	"strconv"
	"strings"

	alignio "github.com/evolbioinfo/goalign/io"
)

// Scanner represents a lexical scanner.
type Scanner struct {
	r *bufio.Reader

	line, prevline         int // current line (for error message)
	position, prevposition int // current position in the line (for error message)
}

// NewScanner returns a new instance of Scanner.
func NewScanner(r io.Reader) *Scanner {
	return &Scanner{r: bufio.NewReader(r), line: 0, prevline: 0, position: 0, prevposition: 0}
}

// reads the next n runes from the bufferred reader.
func (s *Scanner) Read(n int) string {
	var buf bytes.Buffer

	for i := 0; i < 10; i++ {
		ch, _, err := s.r.ReadRune()
		buf.WriteRune(ch)
		if err != nil {
			buf.WriteRune(eof)
			break
		}
		s.prevposition = s.position
		s.position += 1
		if isEndOfLine(ch) {
			s.prevline = s.line
			s.line += 1
			s.position = 0
		}
	}
	return buf.String()
}

// read reads the next rune from the bufferred reader.
// Returns the rune(0) if an error occurs (or io.EOF is returned).
func (s *Scanner) read() rune {
	ch, _, err := s.r.ReadRune()
	s.prevposition = s.position
	s.prevline = s.line
	s.position += 1
	if isEndOfLine(ch) {
		s.line += 1
		s.position = 0
	}

	if err != nil {
		return eof
	}
	return ch
}

// unread places the previously read rune back on the reader.
func (s *Scanner) unread() {
	s.position = s.prevposition
	s.line = s.prevline
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
			}
			alignio.ExitWithMessage(errors.New("\\r without \\n detected"))
		} else {
			return ENDOFLINE, ""
		}
	}

	switch ch {
	case eof:
		return EOF, ""
	}

	s.unread()
	tok, ident := s.scanIdent()

	_, err := strconv.ParseInt(ident, 10, 64)
	if err != nil {
		return tok, ident
	}
	return NUMERIC, ident
}

// scanEndOfLine consumes the current rune and all contiguous \n\r.
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
	switch strings.ToUpper(buf.String()) {
	case "CLUSTAL":
		return CLUSTAL, buf.String()
	case "CLUSTALW":
		return CLUSTAL, buf.String()
	default:
		return IDENTIFIER, buf.String()
	}
}
