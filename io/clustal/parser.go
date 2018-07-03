package clustal

import (
	"bytes"
	"errors"
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/fredericlemoine/goalign/align"
	alignio "github.com/fredericlemoine/goalign/io"
)

// Parser represents a parser.
type Parser struct {
	s      *Scanner
	strict bool
	buf    struct {
		tok Token  // last read token
		lit string // last read literal
		n   int    // buffer size (max=1)
	}
}

// NewParser returns a new instance of Parser.
func NewParser(r io.Reader, strict bool) *Parser {
	return &Parser{s: NewScanner(r), strict: strict}
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
	tok, lit = p.s.Scan(considerWS)

	// Save it to the buffer in case we unscan later.
	p.buf.tok, p.buf.lit = tok, lit
	return
}

func (p *Parser) scanWithEOL() (tok Token, lit string) {
	tok, lit = p.scan()
	if tok != ENDOFLINE {
		return
	}
	prevtok, prevlit := tok, lit
	for ; tok == ENDOFLINE; tok, lit = p.scan() {
		prevtok, prevlit = tok, lit
	}
	p.unscan()
	tok, lit = prevtok, prevlit
	return
}

// unscan pushes the previously read token back onto the buffer.
func (p *Parser) unscan() { p.buf.n = 1 }

// Parse parses a clustal alignment
func (p *Parser) Parse() (align.Alignment, error) {
	var nbseq int64 = 0
	var lenseq int64 = 0
	var al align.Alignment
	var err error
	var nameseqmap map[string]string = make(map[string]string)

	tok, lit := p.scan()
	if tok != CLUSTAL {
		return nil, errors.New("Clustal alignment file must start with 'CLUSTAL'")
	}
	for tok != ENDOFLINE && tok != EOF {
		tok, lit = p.scanWithEOL()
	}

	var name, seq string
	var nbseqs = 0
	var currentnbseqs = 0
	var nblocks = 1
	for tok != EOF {
		// Scan sequence name
		tok, lit = p.scan()
		// Last line of a block
		if tok == WS {
			if currentnbseqs == 0 {
				return nil, errors.New("There should not be a white space here, only at last line of blocks")
			}
			if nbseqs != 0 && currentnbseqs != nbseqs {
				return nil, errors.New(
					fmt.Sprintf("Conservation line: Sequence block nb %d has different number of sequence (%d)",
						nblocks, currentnbseqs))
			}
			for tok != ENDOFLINE && tok != EOF {
				tok, lit = p.scanWithEOL()
			}
			if tok != ENDOFLINE {
				return nil, errors.New("There should be a new line after degree of conservation line")
			}
			tok, lit = p.scanWithEOL()
			if tok == EOF {
				break
			}
			if tok != ENDOFLINE {
				return nil, errors.New("There should be a new line after degree of conservation line")
			}
			tok, lit = p.scan()
			if tok == EOF {
				break
			}
			if nbseqs != 0 && currentnbseqs != nbseqs {
				return nil, errors.New(
					fmt.Sprintf("Sequence block nb %d has different number of sequence (%d)",
						nblocks, currentnbseqs))
			}
			nbseqs = currentnbseqs
			nblocks++
			currentnbseqs = 0
		}

		if tok != IDENTIFIER {
			return nil, errors.New("We expect a sequence identifier here")
		}
		name = lit
		tok, lit = p.scan()
		if tok != WS {
			return nil, errors.New("We expect a whitespace after sequence name")
		}

		tok, lit = p.scan()
		if tok != IDENTIFIER {
			return nil, errors.New("We expect a sequence here")
		}
		seq = lit
		tok, lit = p.scan()
		if tok == NUMERIC {
			// It is a cumulative count of residues
			// skip
			tok, lit = p.scan()
		}
		if tok != ENDOFLINE {
			return nil, errors.New("We expect ENDOFLINE after sequence in a block")
		}
		if tmpseq, ok := nameseqmap[name]; !ok {
			nameseqmap[name] = seq
		} else {
			nameseqmap[name] = tmpseq + seq
		}
	}

	if length(nameseqmap) == 0 {
		return nil, errors.New("No sequences in the alignment")
	}

	var alphabet int = -1
	for k, v := range nameseqmap {
		if al == nil {
			al = align.NewAlign(align.DetectAlphabet(v))
		}
		al.AddSequence(k, v, "")
	}
	return al, nil
}
