package stockholm

import (
	"fmt"
	"io"
	"strings"

	"github.com/evolbioinfo/goalign/align"
)

// Parser represents a parser.
// If ignore is align.IGNORE_NONE: Does not ignore anything
// If ignore is align.IGNORE_NAME: Ignore sequences having the same name (keep the first one whatever their sequence)
// If ignore is align.IGNORE_SEQUENCE: Ignore sequences having the same name and the same sequence
// Otherwise, sets IGNORE_NONE
type Parser struct {
	s               *Scanner
	ignoreidentical int
	alphabet        int // can be align.BOTH, align.AMINOACIDS or align.NUCLEOTIDS
	buf             struct {
		tok Token  // last read token
		lit string // last read literal
		n   int    // buffer size (max=1)
	}
}

// NewParser returns a new instance of Parser.
func NewParser(r io.Reader) *Parser {
	return &Parser{s: NewScanner(r), ignoreidentical: align.IGNORE_NONE, alphabet: align.BOTH}
}

// If sets to true, then will ignore duplicate sequences that have the same name and the same sequence
// Otherwise, it just renames them just as the sequences that have same name and different sequences
func (p *Parser) IgnoreIdentical(ignore int) *Parser {
	p.ignoreidentical = ignore
	return p
}

// alphabet: can be align.BOTH (auto detect alphabet), align.NUCLEOTIDS (considers alignment as nucleotides),
// or align.AMINOACIDS (considers the alignment as aminoacids). If not auto, can return an error if the alignment
// is not compatible with the given alphabet.
// If another value is given, then align.BOTH is considered
func (p *Parser) Alphabet(alphabet int) *Parser {
	p.alphabet = alphabet
	if p.alphabet != align.BOTH &&
		p.alphabet != align.NUCLEOTIDS &&
		p.alphabet != align.AMINOACIDS {
		p.alphabet = align.BOTH
	}
	return p
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

// scanIgnoreWhitespace scans the next non-whitespace token.
func (p *Parser) scanIgnoreWhitespace() (tok Token, lit string) {
	tok, lit = p.scan()
	if tok == WS {
		tok, lit = p.scan()
	}
	return
}

// Parses Stockholm content from the reader
func (p *Parser) Parse() (al align.Alignment, err error) {
	gap := '.'

	// First token should be a "# STOCKHOLM 1.0" token.
	tok, lit := p.scanIgnoreWhitespace()
	if tok != MARKUP {
		err = fmt.Errorf("found %q, expected # STOCKHOLM 1.0", lit)
		return
	}
	tok, lit = p.scanIgnoreWhitespace()
	if tok != STOCKHOLM {
		err = fmt.Errorf("found %q, expected # STOCKHOLM 1.0", lit)
		return
	}
	_, lit = p.scanIgnoreWhitespace()
	if lit != "1.0" {
		err = fmt.Errorf("found %q, expected # STOCKHOLM 1.0", lit)
		return
	}

	al = align.NewAlign(align.UNKNOWN)
	al.IgnoreIdentical(p.ignoreidentical)

	// Now we can parse the remaining of the file
	for {
		tok, lit := p.scanIgnoreWhitespace()
		if tok == ILLEGAL {
			err = fmt.Errorf("found illegal token %q", lit)
			return
		}
		if tok == EOF {
			break
		}
		if tok == ENDOFLINE {
			continue
		}

		if tok == MARKUP {
			for tok != ENDOFLINE {
				tok, _ = p.scanIgnoreWhitespace()
			}
			continue
		}

		if tok == END {
			break
		}

		if tok == IDENT || tok == NUMERIC {
			name := lit
			tok, lit = p.scanIgnoreWhitespace()
			if tok != IDENT {
				err = fmt.Errorf("found illegal sequence %q", lit)
				return
			}
			sequence := lit
			sequence = strings.Replace(sequence, string(gap), string(align.GAP), -1)
			if err = al.AddSequence(name, sequence, ""); err != nil {
				return
			}
		}
	}

	if al.Length() == 0 {
		err = fmt.Errorf("no sequence in this Stockholm file")
		return
	}

	// If the alphabet given on the command line is BOTH,
	// then we take the alphabet given in the stockholm file
	if p.alphabet == align.BOTH {
		al.AutoAlphabet()
	} else {
		err = al.SetAlphabet(p.alphabet)
		return
	}
	return
}
