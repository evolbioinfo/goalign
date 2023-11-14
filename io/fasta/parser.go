package fasta

import (
	"bytes"
	"errors"
	"io"
	"regexp"
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
func (p *Parser) Parse() (al align.Alignment, err error) {
	al = align.NewAlign(align.UNKNOWN)
	al.IgnoreIdentical(p.ignoreidentical)
	err = p.parseGeneric(al)
	return
}

// ParseUnalign parses an unaligned FASTA file
func (p *Parser) ParseUnalign() (sb align.SeqBag, err error) {
	sb = align.NewSeqBag(align.UNKNOWN)
	sb.IgnoreIdentical(p.ignoreidentical)
	err = p.parseGeneric(sb)
	return
}

func (p *Parser) parseGeneric(sb align.SeqBag) (err error) {
	// The first token should be a ">"
	tok, lit := p.scanIgnoreEndOfLine()
	if tok != STARTIDENT {
		err = errors.New("fasta file should start with a > ")
		return
	}
	p.unscan()
	var curseq bytes.Buffer
	curname := ""
	firstSpaces := regexp.MustCompile("^( +)")
	for tok != EOF {
		tok, lit = p.scanIgnoreEndOfLine()
		switch tok {
		case STARTIDENT:
			tok, lit = p.scan()
			if tok != IDENTIFIER {
				err = errors.New("> should be followed by a sequence identifier")
				return
			}
			if curseq.Len() > 0 {
				if err = sb.AddSequence(curname, curseq.String(), ""); err != nil {
					return
				}
				curseq.Reset()
			} else if curname != "" {
				err = errors.New("A Fasta entry has a name but no sequence (" + curname + ")")
				return
			}
			curname = firstSpaces.ReplaceAllString(lit, "")
		case IDENTIFIER:
			curseq.WriteString(strings.Replace(lit, " ", "", -1))
		case EOF:
			if curseq.Len() > 0 {
				if err = sb.AddSequence(curname, curseq.String(), ""); err != nil {
					return
				}
			}
		}
	}

	if p.alphabet == align.BOTH {
		sb.AutoAlphabet()
	} else {
		if err = sb.SetAlphabet(p.alphabet); err != nil {
			return
		}
	}
	return
}
