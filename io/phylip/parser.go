package phylip

import (
	"bytes"
	"errors"
	"fmt"
	"github.com/fredericlemoine/goalign/align"
	alignio "github.com/fredericlemoine/goalign/io"
	"io"
	"strconv"
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

// Parse parses a phylip alignment
func (p *Parser) Parse() (align.Alignment, error) {
	var nbseq int64 = 0
	var lenseq int64 = 0
	var al align.Alignment = nil
	var err error

	// We skip all WS and EOL at the beginning
	tok, lit := p.scanWithEOL()
	for tok == WS || tok == ENDOFLINE {
		tok, lit = p.scanWithEOL()
	}
	if tok == EOF {
		return nil, nil
	}

	// The first token different from WS and EOL should be a Number
	if tok != NUMERIC {
		return nil, errors.New("Phylip file must begin with the number of sequences : " + fmt.Sprintf("%d", tok))
	} else {
		nbseq, err = strconv.ParseInt(lit, 10, 64)
		if err != nil {
			return nil, errors.New("The numeric is not parsable: " + lit)
		}
	}

	names := make([]string, nbseq)
	seqs := make([]*bytes.Buffer, nbseq)

	tok, lit = p.scan()
	if tok != WS {
		return nil, errors.New("There should be a whitespace between number of sequences and length")
	}

	// The second token, after a WS, should be a Number
	tok, lit = p.scan()
	if tok != NUMERIC {
		return nil, errors.New("Phylip file must begin with the number of sequences and their length")
	} else {
		lenseq, err = strconv.ParseInt(lit, 10, 64)
		if err != nil {
			return nil, errors.New("The numeric is not parsable: " + lit)
		}
	}

	// Then a \n
	tok, lit = p.scan()
	if tok != ENDOFLINE {
		return nil, errors.New("Bad Phylip format, \n missing after header")
	}

	// Then names of the sequences and sequences
	for i := 0; i < int(nbseq); i++ {
		tok, lit = p.scan()
		if tok != IDENTIFIER && tok != NUMERIC {
			return nil, errors.New("Bad Phylip format, we should have an sequence identifier after the header : " + lit)
		}
		names[i] = lit
		tok, lit = p.scan()
		seqs[i] = new(bytes.Buffer)
		for tok != ENDOFLINE {
			switch tok {
			case IDENTIFIER:
				seqs[i].WriteString(lit)
			case WS:
			default:
				return nil, errors.New("Bad Phylip format, Unexpected character :" + lit)
			}
			tok, lit = p.scan()
		}
	}

	tok, lit = p.scanWithEOL()

	if tok == ENDOFLINE {
		tok, lit = p.scan()
		p.unscan()
	} else if int(lenseq) != seqs[0].Len() {
		alignio.ExitWithMessage(errors.New("Bad Phylip Format : Should have a blank line here"))
	}
	//  else if tok != EOF {
	// 	alignio.ExitWithMessage(errors.New("Bad Phylip Format : Should not have a character here, all sequences have been red"))
	// }

	// Then other blocks with only sequences
	b := 0
	for tok != EOF && int(lenseq) != seqs[0].Len() {
		for i := 0; i < int(nbseq); i++ {
			tok, lit := p.scan()
			if tok == WS {
				tok, lit = p.scan()
			}

			if tok != IDENTIFIER {
				fmt.Println("Block: ")
				fmt.Println(b)
				return nil, errors.New("Bad Phylip format, we should have a sequence block here")
			}
			for tok != ENDOFLINE {
				switch tok {
				case IDENTIFIER:
					seqs[i].WriteString(lit)
				case WS:
				default:
					return nil, errors.New("Bad Phylip format, Unexpected character :" + lit)
				}
				tok, lit = p.scan()
			}
		}
		tok, lit = p.scanWithEOL()
		if tok == ENDOFLINE {
			tok, lit = p.scan()
			p.unscan()
		} else if int(lenseq) != seqs[0].Len() {
			alignio.ExitWithMessage(errors.New("Bad Phylip Format : Should have a blank line here"))
		}
		b++
	}
	tok, lit = p.scanWithEOL()
	p.unscan()

	for i, name := range names {
		seq := seqs[i].String()
		if int(lenseq) != len(seq) {
			alignio.ExitWithMessage(errors.New("Bad Phylip format : Length of sequences does not correspond to header"))
		}
		if al == nil {
			al = align.NewAlign(align.DetectAlphabet(seq))
		}
		al.AddSequence(name, seq, "")
	}
	return al, nil
}

/*
Parses mutliple phylip alignments and give them to the channel
At the end (or in case of error), it closes the channel
*/
func (p *Parser) ParseMultiple(aligns chan<- align.Alignment) error {
	al, err := p.Parse()
	if err != nil {
		close(aligns)
		return err
	}
	for err == nil && al != nil {
		aligns <- al
		al, err = p.Parse()
		if err != nil {
			close(aligns)
			return err
		}
	}
	close(aligns)
	return nil
}
