package partition

import (
	"errors"
	"fmt"
	"io"
	"strconv"

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

// unscan pushes the previously read token back onto the buffer.
func (p *Parser) unscan() { p.buf.n = 1 }

// Parse parses a partition file
func (p *Parser) Parse(alignmentLength int) (partitionSet *align.PartitionSet, err error) {
	partitionSet = align.NewPartitionSet(alignmentLength)
	err = p.parse(partitionSet)
	return
}

func (p *Parser) parse(ps *align.PartitionSet) (err error) {
	// The first token should be a ">"
	tok, lit := p.scan()
	if tok != IDENTIFIER {
		err = errors.New("Partition should start with a Model name")
		return
	}
	p.unscan()

	var start, end, modulo int64
	var modeleName, partitionName string

	for tok != EOF {
		tok, lit = p.scan()
		switch tok {
		case IDENTIFIER:
			modeleName = lit
			tok, lit = p.scan()
			if tok != SEPARATOR {
				err = fmt.Errorf("Modele name should be followed by ',' : %s", modeleName)
				return
			}

			tok, lit = p.scan()
			partitionName = lit
			if tok != IDENTIFIER {
				err = fmt.Errorf("Modele name should be followed by ',' then the name of the partition: %s", lit)
				return
			}

			tok, lit = p.scan()
			if tok != EQUAL {
				err = fmt.Errorf("Partition name should be followed by '=' : [%s|%s]", partitionName, lit)
				return
			}
			// Parse intervals
			tok, lit = p.scan()
			for tok != ENDOFLINE && tok != EOF {
				if tok != DECIMAL {
					err = fmt.Errorf("Interval definition should start with a number : [%s]", lit)
					return
				}
				start, _ = strconv.ParseInt(lit, 10, 64)
				end = start
				tok, lit = p.scan()
				if tok == RANGE {
					tok, lit = p.scan()
					if tok != DECIMAL {
						err = fmt.Errorf("Interval definition '-' should be followed by an integer value : [%s]", lit)
						return
					}
					end, _ = strconv.ParseInt(lit, 10, 64)
					tok, lit = p.scan()
				}
				if tok == MODULO {
					tok, lit = p.scan()
					if tok != DECIMAL {
						err = fmt.Errorf("there should be an integer value after '/': [%s]", lit)
						return
					}
					modulo, _ = strconv.ParseInt(lit, 10, 64)
					tok, lit = p.scan()
				} else {
					modulo = 1
				}

				if err = ps.AddRange(partitionName, modeleName, int(start)-1, int(end)-1, int(modulo)); err != nil {
					return
				}

				if tok == SEPARATOR {
					tok, lit = p.scan()
				} else if tok != ENDOFLINE && tok != EOF {
					err = fmt.Errorf("there should be a separator (or EOL or EOF) after interval definition : [%s]", lit)
					return
				}
			}
		}
	}
	return
}
