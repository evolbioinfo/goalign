package clustal

import (
	"errors"
	"fmt"
	"io"

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
func (p *Parser) Parse() (al align.Alignment, err error) {
	var nbseq int = 0
	var names []string = make([]string, 0)
	var seqs []string = make([]string, 0)
	tok, lit := p.scan()
	if tok != CLUSTAL {
		return nil, errors.New("clustal alignment file must start with 'CLUSTAL'")
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
				return nil, errors.New("there should not be a white space here, only at last line of blocks")
			}
			if nbseq != 0 && currentnbseqs != nbseq {
				err = fmt.Errorf("conservation line: Sequence block nb %d has different number of sequence (%d)",
					nblocks, currentnbseqs)
				return
			}
			for tok != ENDOFLINE && tok != EOF {
				tok, lit = p.scan()
			}
			if tok != ENDOFLINE {
				err = fmt.Errorf("there should be a new line after degree of conservation line, line  %d, position %d, have \"%s\", block %d, sequence %d", p.s.line+1, p.s.position+1, lit, nblocks, currentnbseqs)
				return
			}
			tok, lit = p.scanWithEOL()
			if tok == EOF {
				break
			}
			if tok != ENDOFLINE {
				err = fmt.Errorf("expecting a new line after sequence conservation line,  line %d, position %d, have \"%s\", block %d, sequence %d", p.s.line+1, p.s.position+1, lit, nblocks, currentnbseqs)
				return
			}
			tok, lit = p.scan()
			if tok == EOF {
				break
			}
			if nbseq != 0 && currentnbseqs != nbseq {
				err = fmt.Errorf("sequence block has different number of sequence, line %d, position %d, have \"%s\", block %d, sequence %d", p.s.line+1, p.s.position+1, lit, nblocks, currentnbseqs)
				return
			}
			nbseq = currentnbseqs
			nblocks++
			currentnbseqs = 0
		}

		if tok != IDENTIFIER && tok != NUMERIC {
			err = fmt.Errorf("expecting a sequence identifier (or a sequence conservation line) line %d, position %d, have \"%s\", block %d, sequence %d", p.s.line+1, p.s.position+1, lit, nblocks, currentnbseqs)
			return
		}
		name = lit
		tok, lit = p.scan()
		if tok != WS {
			err = fmt.Errorf("expecting a whitespace after sequence name,  line %d, position %d, have \"%s\", block %d, sequence %d", p.s.line+1, p.s.position+1, lit, nblocks, currentnbseqs)
			return
		}

		tok, lit = p.scan()
		if tok != IDENTIFIER {
			err = fmt.Errorf("expecting a sequence here, line %d, position %d, have \"%s\", block %d, sequence %d", p.s.line+1, p.s.position+1, lit, nblocks, currentnbseqs)
			return
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
				err = fmt.Errorf("expecting a current length after sequence + whitespace, line %d, position %d, have \"%s\", block %d, sequence %d", p.s.line+1, p.s.position+1, lit, nblocks, currentnbseqs)
				return
			}
		}
		if tok != ENDOFLINE {
			err = fmt.Errorf("expecting ENDOFLINE after sequence in a block, line %d, position %d, have \"%s\", block %d, sequence %d", p.s.line+1, p.s.position+1, lit, nblocks, currentnbseqs)
			return
		}

		if nblocks == 0 {
			names = append(names, name)
			seqs = append(seqs, seq)
		} else {
			if names[currentnbseqs] != name {
				err = fmt.Errorf("sequence name in the current block (%s) does not correspond to sequence name in first block (%s),  line %d, position %d, have \"%s\", block %d, sequence %d", names[currentnbseqs], name, p.s.line+1, p.s.position+1, lit, nblocks, currentnbseqs)
				return
			}
			seqs[currentnbseqs] = seqs[currentnbseqs] + seq
		}
		currentnbseqs++
	}

	if len(names) == 0 {
		err = errors.New("no sequences in the alignment")
		return
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

	if p.alphabet == align.BOTH {
		al.AutoAlphabet()
	} else {
		if err = al.SetAlphabet(p.alphabet); err != nil {
			return
		}
	}

	return
}
