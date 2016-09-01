package fasta

import (
	"bytes"
	"errors"
	"github.com/fredericlemoine/goalign/align"
	"io"
)

// Parser represents a parser.
type Parser struct {
	s   *Scanner
	buf struct {
		tok Token  // last read token
		lit string // last read literal
		n   int    // buffer size (max=1)
	}
}

// NewParser returns a new instance of Parser.
func NewParser(r io.Reader) *Parser {
	return &Parser{s: NewScanner(r)}
}

// scan returns the next token from the underlying scanner.
// If a token has been unscanned then read that instead.
func (p *Parser) scan() (tok Token, lit string) {
	// If we have a token on the buffer, then return it.
	if p.buf.n != 0 {
		p.buf.n = 0
		return p.buf.tok, p.buf.lit
	}

	// Otherwise read the next token from the scanner.
	tok, lit = p.s.Scan()

	// Save it to the buffer in case we unscan later.
	p.buf.tok, p.buf.lit = tok, lit
	return
}

// unscan pushes the previously read token back onto the buffer.
func (p *Parser) unscan() { p.buf.n = 1 }

// scanIgnoreEndOfLine scans the next EOL token.
func (p *Parser) scanIgnoreEndOfLine() (tok Token, lit string) {
	tok, lit = p.scan()
	if tok == ENDOFLINE {
		tok, lit = p.scan()
	}
	return
}

// Parse parses a FASTA alignment
func (p *Parser) Parse() (align.Alignment, error) {
	// The first token should be a ">"
	tok, lit := p.scanIgnoreEndOfLine()
	if tok != STARTIDENT {
		return nil, errors.New("Fasta file should start with a > ")
	}
	p.unscan()
	var al align.Alignment = nil
	var curseq bytes.Buffer
	curname := ""

	for tok != EOF {
		tok, lit = p.scanIgnoreEndOfLine()

		switch tok {
		case STARTIDENT:
			tok, lit = p.scan()
			if tok != IDENTIFIER {
				return nil, errors.New("> should be followed by a sequence identifier")
			}
			if curseq.Len() > 0 {
				s := curseq.String()
				if al == nil {
					al = align.NewAlign(align.DetectAlphabet(s))
				}
				err := al.AddSequence(curname, curseq.String(), "")
				if err != nil {
					return nil, err
				}
				curseq.Reset()
			} else if curname != "" {
				return nil, errors.New("A Fasta entry has a name but no sequence (" + curname + ")")
			}
			curname = lit
		case IDENTIFIER:
			curseq.WriteString(lit)
		case EOF:
			if curseq.Len() > 0 {
				s := curseq.String()
				if al == nil {
					al = align.NewAlign(align.DetectAlphabet(s))
				}
				err := al.AddSequence(curname, s, "")
				if err != nil {
					return nil, err
				}
			}
		}
	}
	return al, nil
}
