package phylip

import (
	"bytes"
	"errors"
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/evolbioinfo/goalign/align"
)

// Parser represents a parser.
type Parser struct {
	s               *Scanner
	strict          bool
	ignoreidentical int
	alphabet        int // can be align.BOTH, align.AMINOACIDS or align.NUCLEOTIDS
	buf             struct {
		tok Token  // last read token
		lit string // last read literal
		n   int    // buffer size (max=1)
	}
}

// NewParser returns a new instance of Parser.
func NewParser(r io.Reader, strict bool) *Parser {
	return &Parser{s: NewScanner(r), strict: strict, ignoreidentical: align.IGNORE_NONE, alphabet: align.BOTH}
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

// Parse parses a phylip alignment
func (p *Parser) Parse() (al align.Alignment, err error) {
	var nbseq int64 = 0
	var lenseq int64 = 0

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
		if nbseq == 0 {
			err = fmt.Errorf("no sequences in the alignment")
			return
		}
		if nbseq < 0 {
			err = fmt.Errorf("wrong number of sequences in the alignment: %d", nbseq)
			return
		}
	}

	names := make([]string, nbseq)
	seqs := make([]*bytes.Buffer, nbseq)

	tok, lit = p.scan()
	if tok != WS {
		err = fmt.Errorf("there should be a whitespace between number of sequences and length")
		return
	}

	// The second token, after a WS, should be a Number
	tok, lit = p.scan()
	if tok != NUMERIC {
		err = fmt.Errorf("phylip file must begin with the number of sequences and their length")
		return
	} else {
		lenseq, err = strconv.ParseInt(lit, 10, 64)
		if err != nil {
			err = fmt.Errorf("the numeric is not parsable: %s", lit)
			return
		}
		if lenseq == 0 {
			err = fmt.Errorf("0 length sequences defined in the header")
			return
		}
		if nbseq < 0 {
			err = fmt.Errorf("wrong sequence length in the header: %d", lenseq)
			return
		}
	}

	// Then a \n
	tok, lit = p.scan()
	if tok != ENDOFLINE {
		err = fmt.Errorf("bad Phylip format, \\n missing after header")
		return
	}

	// Then names of the sequences and sequences
	for i := 0; i < int(nbseq); i++ {
		if p.strict {
			// if strict
			name := p.s.Read(10)
			if []rune(name)[len(name)-1] == eof {
				err = fmt.Errorf("bad Phylip format, less sequences in the file than indicated in the header : %d vs. %d", nbseq, i)
				return
			}
			if name == "" {
				err = fmt.Errorf("bad Phylip format, we should have an sequence identifier after the header : %s", lit)
				return
			}
			// We remove spaces from names...
			name = strings.Replace(name, " ", "", -1)
			names[i] = name
		} else {
			tok, lit = p.scan()
			if tok == EOF {
				err = fmt.Errorf("bad Phylip format, less sequences in the file than indicated in the header : %d vs. %d", nbseq, i)
				return
			}
			if tok != IDENTIFIER && tok != NUMERIC {
				err = fmt.Errorf("bad Phylip format, we should have an sequence identifier after the header : %s", lit)
				return
			}
			names[i] = lit
		}

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
		err = fmt.Errorf("bad Phylip Format : Should have a blank line here")
		return
	}
	// All sequences are completely parsed
	// If there are several alignments in the file, we should
	// unscan the last token
	if int(lenseq) == seqs[0].Len() {
		p.unscan()
	}

	// Then other blocks with only sequences
	// If the sequences have not been already complety parsed
	b := 0
	for tok != EOF && int(lenseq) != seqs[0].Len() {
		for i := 0; i < int(nbseq); i++ {
			tok, lit := p.scan()
			if tok == WS {
				tok, lit = p.scan()
			}

			if tok != IDENTIFIER {
				// fmt.Println("Block: ")
				// fmt.Println(b)
				err = fmt.Errorf("bad Phylip format, we should have a sequence block here")
				return
			}
			for tok != ENDOFLINE {
				switch tok {
				case IDENTIFIER:
					seqs[i].WriteString(lit)
				case WS:
				default:
					err = fmt.Errorf("bad Phylip format, unexpected character : %s", lit)
					return
				}
				tok, lit = p.scan()
			}
		}
		tok, lit = p.scanWithEOL()
		if tok == ENDOFLINE {
			tok, lit = p.scan()
			p.unscan()
		} else if int(lenseq) != seqs[0].Len() {
			err = fmt.Errorf("bad Phylip Format : Should have a blank line here")
			return
		}
		b++
	}
	//tok, lit = p.scanWithEOL()
	//p.unscan()

	al = align.NewAlign(align.UNKNOWN)
	al.IgnoreIdentical(p.ignoreidentical)
	for i, name := range names {
		seq := seqs[i].String()
		if int(lenseq) != len(seq) {
			err = fmt.Errorf("bad Phylip format : Length of sequences does not correspond to header")
			return
		}
		al.AddSequence(name, seq, "")
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

/*
Parses mutliple phylip alignments and give them to the channel
At the end (or in case of error), it closes the channel.
It is better to call this functioni from a go routine
since it can block the buffered channel.
*/
func (p *Parser) ParseMultiple(aligns *align.AlignChannel) {
	var al align.Alignment
	var err error
	al, err = p.Parse()
	for err == nil && al != nil {
		aligns.Achan <- al
		al, err = p.Parse()
	}
	aligns.Err = err

	close(aligns.Achan)
}
