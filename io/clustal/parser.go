package clustal

import (
	"errors"
	"fmt"
	"io"

	"github.com/evolbioinfo/goalign/align"
)

// Parser represents a parser.
type Parser struct {
	s               *Scanner
	ignoreidentical bool
	buf             struct {
		tok Token  // last read token
		lit string // last read literal
		n   int    // buffer size (max=1)
	}
}

// NewParser returns a new instance of Parser.
func NewParser(r io.Reader) *Parser {
	return &Parser{s: NewScanner(r), ignoreidentical: false}
}

// If sets to true, then will ignore duplicate sequences that have the same name and the same sequence
// Otherwise, it just renames them just as the sequences that have same name and different sequences
func (p *Parser) IgnoreIdentical(ignore bool) {
	p.ignoreidentical = ignore
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
	var names []string = make([]string, 0)
	var seqs []string = make([]string, 0)
	tok, lit := p.scan()
	if tok != CLUSTAL {
		return nil, errors.New("Clustal alignment file must start with 'CLUSTAL'")
	}
	for tok != ENDOFLINE && tok != EOF {
		tok, lit = p.scanWithEOL()
	}

	var name, seq string
	var currentnbseqs int = 0
	var nblocks = 0
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

		if nblocks == 0 {
			names = append(names, name)
			seqs = append(seqs, seq)
		} else {
			if names[currentnbseqs] != name {
				return nil, fmt.Errorf("Name at block %d line %d does not correspond to name in first block (%s vs. %s)", nblocks, currentnbseqs, names[currentnbseqs], name)
			}
			seqs[currentnbseqs] = seqs[currentnbseqs] + seq
		}
		currentnbseqs++
	}

	if len(names) == 0 {
		return nil, errors.New("No sequences in the alignment")
	}

	for i, n := range names {
		s := seqs[i]
		if al == nil {
			al = align.NewAlign(align.DetectAlphabet(s))
			al.IgnoreIdentical(p.ignoreidentical)
		}
		if err := al.AddSequence(n, s, ""); err != nil {
			return nil, err
		}
	}

	return al, nil
}
