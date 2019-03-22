package clustal

import (
	"errors"
	"fmt"
	"io"

	"github.com/evolbioinfo/goalign/align"
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
	var nbseq int = 0
	var al align.Alignment
	var nameseqmap map[string]string = make(map[string]string)
	var seqnames []string = make([]string, 0)
	tok, lit := p.scan()
	if tok != CLUSTAL {
		return nil, errors.New("Clustal alignment file must start with 'CLUSTAL'")
	}
	for tok != ENDOFLINE && tok != EOF {
		tok, lit = p.scanWithEOL()
	}

	var name, seq string
	var currentnbseqs int = 0
	var nblocks = 1
	for tok != EOF {
		// Scan sequence name
		tok, lit = p.scan()
		// Last line of a block
		if tok == WS {
			if currentnbseqs == 0 {
				return nil, errors.New("There should not be a white space here, only at last line of blocks")
			}
			if nbseq != 0 && currentnbseqs != nbseq {
				return nil, errors.New(
					fmt.Sprintf("Conservation line: Sequence block nb %d has different number of sequence (%d)",
						nblocks, currentnbseqs))
			}
			for tok != ENDOFLINE && tok != EOF {
				tok, lit = p.scan()
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
			if nbseq != 0 && currentnbseqs != nbseq {
				return nil, errors.New(
					fmt.Sprintf("Sequence block nb %d has different number of sequence (%d)",
						nblocks, currentnbseqs))
			}
			nbseq = currentnbseqs
			nblocks++
			currentnbseqs = 0
		}

		if tok != IDENTIFIER && tok != NUMERIC {
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
		currentnbseqs++
		tok, lit = p.scan()
		if tok == WS {
			tok, lit = p.scan()
			if tok == NUMERIC {
				// It is a cumulative count of residues
				// skip
				tok, lit = p.scan()
			} else {
				return nil, errors.New("We expect a current length after sequence + whitespace")
			}
		}
		if tok != ENDOFLINE {
			return nil, errors.New("We expect ENDOFLINE after sequence in a block")
		}
		if tmpseq, ok := nameseqmap[name]; !ok {
			if nblocks > 1 {
				return nil, errors.New(
					fmt.Sprintf(
						"There should not be a new sequence name (%s) found in block nb %d",
						name,
						nblocks,
					))
			}
			nameseqmap[name] = seq
			seqnames = append(seqnames, name)
		} else {
			nameseqmap[name] = tmpseq + seq
		}
	}

	if len(nameseqmap) == 0 {
		return nil, errors.New("No sequences in the alignment")
	}

	for _, n := range seqnames {
		if s, ok := nameseqmap[n]; !ok {
			return nil, errors.New("Error in clustal alignment parsing: no sequence with name " + n)
		} else {
			if al == nil {
				al = align.NewAlign(align.DetectAlphabet(s))
			}
			if err := al.AddSequence(n, s, ""); err != nil {
				return nil, err
			}
		}
	}

	return al, nil
}
