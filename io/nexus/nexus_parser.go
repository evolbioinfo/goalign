package nexus

import (
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/evolbioinfo/goalign/align"
	aio "github.com/evolbioinfo/goalign/io"
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

// unscan pushes the previously read token back onto the buffer.
func (p *Parser) unscan() { p.buf.n = 1 }

// scanIgnoreWhitespace scans the next non-whitespace token.
func (p *Parser) scanIgnoreWhitespace() (tok Token, lit string) {
	tok, lit = p.scan()
	if tok == WS {
		tok, lit = p.scan()
	}
	return
}

// Parses Nexus content from the reader
func (p *Parser) Parse() (al align.Alignment, err error) {
	var nchar, ntax, taxantax int64
	datatype := "dna"
	missing := '*'
	matchchar := '.'
	gap := '-'
	var taxlabels map[string]bool = nil
	var names []string
	var sequences map[string]string

	// First token should be a "NEXUS" token.
	tok, lit := p.scanIgnoreWhitespace()
	if tok != NEXUS {
		err = fmt.Errorf("found %q, expected #NEXUS", lit)
		return
	}

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

		if tok == OPENBRACK {
			if tok, lit, err = p.consumeComment(tok, lit); err != nil {
				return
			}
			tok, lit = p.scanIgnoreWhitespace()
		}

		// Beginning of a block
		if tok == BEGIN {
			// Next token should be the name of the block
			tok2, lit2 := p.scanIgnoreWhitespace()
			// Then a ;
			tok3, lit3 := p.scanIgnoreWhitespace()
			if tok3 != ENDOFCOMMAND {
				err = fmt.Errorf("found %q, expected ;", lit3)
				return
			}
			switch tok2 {
			case TAXA:
				// TAXA BLOCK
				taxantax, taxlabels, err = p.parseTaxa()
			case TREES:
				// If an unsupported block is seen, we just skip it
				aio.PrintMessage("TREE block is not supported in goalign nexus parser, see gotree nexus parser, skipping")
				err = p.parseUnsupportedBlock()
			case DATA:
				// DATA/CHARACTERS BLOCK
				names, sequences, nchar, ntax, datatype, missing, gap, matchchar, err = p.parseData()
			default:
				// If an unsupported block is seen, we just skip it
				aio.PrintMessage(fmt.Sprintf("Unsupported block %q, skipping", lit2))
				err = p.parseUnsupportedBlock()
			}

			if err != nil {
				return
			}
		}
	}

	if int(taxantax) != -1 && int(taxantax) != len(taxlabels) {
		err = fmt.Errorf("Number of defined taxa in TAXLABELS/DIMENSIONS (%d) is different from length of taxa list (%d)", taxantax, len(taxlabels))
		return
	}

	// if gap != '-' || missing != '*' {
	// 	err = fmt.Errorf("We only accept - gaps (have %c) && * missing (have %c) so far", gap, missing)
	// 	return
	// }

	if sequences == nil || len(sequences) == 0 {
		err = fmt.Errorf("No sequence in this Nexus file")
		return
	}

	// We initialize alignment structure using goalign structure
	if names != nil && sequences != nil {
		al = align.NewAlign(align.AlphabetFromString(datatype))
		al.IgnoreIdentical(p.ignoreidentical)
		if al.Alphabet() == align.UNKNOWN {
			err = fmt.Errorf("Unknown datatype: %q", datatype)
			return
		}
		if len(names) != int(ntax) && ntax != -1 {
			err = fmt.Errorf("Number of taxa in alignment (%d)  does not correspond to definition %d", len(names), ntax)
			return
		}
		for i, name := range names {
			seq, _ := sequences[name]
			if len(seq) != int(nchar) && nchar != -1 {
				err = fmt.Errorf("Number of character in sequence #%d (%d) does not correspond to definition %d", i, len(seq), nchar)
				return
			}
			seq = strings.Replace(seq, string(gap), string(align.GAP), -1)
			seq = strings.Replace(seq, string(missing), string(align.OTHER), -1)
			seq = strings.Replace(seq, string(matchchar), string(align.POINT), -1)
			if err = al.AddSequence(name, seq, ""); err != nil {
				return
			}
		}
		// We check that tax labels are the same as alignment sequence names
		if taxlabels != nil {
			al.Iterate(func(name string, sequence string) bool {
				if _, ok := taxlabels[name]; !ok {
					err = fmt.Errorf("Sequence name %s in the alignment is not defined in the TAXLABELS block", name)
				}
				return false
			})
			if err != nil {
				return nil, err
			}
			if al.NbSequences() != len(taxlabels) {
				err = fmt.Errorf("Some taxa names defined in TAXLABELS are not present in the alignment")
				return
			}
		}
		al.ReplaceMatchChars()
	}
	return
}

// Parse taxa block
func (p *Parser) parseTaxa() (int64, map[string]bool, error) {
	taxlabels := make(map[string]bool)
	var err error
	stoptaxa := false
	var ntax int64 = -1
	for !stoptaxa {
		tok, lit := p.scanIgnoreWhitespace()
		switch tok {
		case ENDOFLINE:
			continue
		case ILLEGAL:
			err = fmt.Errorf("found illegal token %q", lit)
			stoptaxa = true
		case EOF:
			err = fmt.Errorf("End of file within a TAXA block (no END;)")
			stoptaxa = true
		case END:
			tok2, _ := p.scanIgnoreWhitespace()
			if tok2 != ENDOFCOMMAND {
				err = fmt.Errorf("End token without ;")
			}
			stoptaxa = true
		case DIMENSIONS:
			// Dimensions of the data: ntax
			stopdimensions := false
			for !stopdimensions {
				tok2, lit2 := p.scanIgnoreWhitespace()
				switch tok2 {
				case ENDOFCOMMAND:
					stopdimensions = true
				case NTAX:
					tok3, lit3 := p.scanIgnoreWhitespace()
					if tok3 != EQUAL {
						err = fmt.Errorf("Expecting '=' after NTAX, got %q", lit3)
						stopdimensions = true
					}
					tok4, lit4 := p.scanIgnoreWhitespace()
					if tok4 != NUMERIC {
						err = fmt.Errorf("Expecting Integer value after 'NTAX=', got %q", lit4)
						stopdimensions = true
					}
					ntax, err = strconv.ParseInt(lit4, 10, 64)
					if err != nil {
						stopdimensions = true
					}
				default:
					if err = p.parseUnsupportedKey(lit2); err != nil {
						stopdimensions = true
					}
					aio.PrintMessage(fmt.Sprintf("Unsupported key %q in %q command, skipping", lit2, lit))
				}
				if err != nil {
					stoptaxa = true
				}
			}
		case TAXLABELS:
			stoplabels := false
			for !stoplabels {
				tok2, lit2 := p.scanIgnoreWhitespace()
				switch tok2 {
				case ENDOFCOMMAND:
					stoplabels = true
				case IDENT:
					taxlabels[lit2] = true
				case ENDOFLINE:
					continue
				default:
					err = fmt.Errorf("Unknown token %q (%v) in taxlabel list", lit2, tok2)
					stoplabels = true
				}
			}
			if err != nil {
				stoptaxa = true
			}
		case OPENBRACK:
			if tok, lit, err = p.consumeComment(tok, lit); err != nil {
				stoptaxa = true
			}

		default:
			err = p.parseUnsupportedCommand()
			aio.PrintMessage(fmt.Sprintf("Unsupported command %q in block TAXA, skipping", lit))
			if err != nil {
				stoptaxa = true
			}
		}
	}
	return ntax, taxlabels, err
}

// DATA / Characters BLOCK
func (p *Parser) parseData() (names []string, sequences map[string]string, nchar, ntax int64, datatype string, missing, gap, matchchar rune, err error) {
	datatype = "dna"
	missing = '*'
	matchchar = '.'
	gap = '-'
	stopdata := false
	sequences = make(map[string]string)
	names = make([]string, 0)
	nchar = -1
	ntax = -1
	for !stopdata {
		tok, lit := p.scanIgnoreWhitespace()
		switch tok {
		case ENDOFLINE:
			break
		case ILLEGAL:
			err = fmt.Errorf("found illegal token %q", lit)
			stopdata = true
		case EOF:
			err = fmt.Errorf("End of file within a TAXA block (no END;)")
			stopdata = true
		case END:
			tok2, _ := p.scanIgnoreWhitespace()
			if tok2 != ENDOFCOMMAND {
				err = fmt.Errorf("End token without ;")
			}
			stopdata = true
		case DIMENSIONS:
			// Dimensions of the data: nchar , ntax
			stopdimensions := false
			for !stopdimensions {
				tok2, lit2 := p.scanIgnoreWhitespace()
				switch tok2 {
				case ENDOFCOMMAND:
					stopdimensions = true
				case NTAX:
					tok3, lit3 := p.scanIgnoreWhitespace()
					if tok3 != EQUAL {
						err = fmt.Errorf("Expecting '=' after NTAX, got %q", lit3)
						stopdimensions = true
					}
					tok4, lit4 := p.scanIgnoreWhitespace()
					if tok4 != NUMERIC {
						err = fmt.Errorf("Expecting Integer value after 'NTAX=', got %q", lit4)
						stopdimensions = true
					}
					ntax, err = strconv.ParseInt(lit4, 10, 64)
					if err != nil {
						stopdimensions = true
					}
				case NCHAR:
					tok3, lit3 := p.scanIgnoreWhitespace()
					if tok3 != EQUAL {
						err = fmt.Errorf("Expecting '=' after NTAX, got %q", lit3)
						stopdimensions = true
					}
					tok4, lit4 := p.scanIgnoreWhitespace()
					if tok4 != NUMERIC {
						err = fmt.Errorf("Expecting Integer value after 'NTAX=', got %q", lit4)
						stopdimensions = true
					}
					nchar, err = strconv.ParseInt(lit4, 10, 64)
					if err != nil {
						stopdimensions = true
					}
				default:
					if err = p.parseUnsupportedKey(lit2); err != nil {
						stopdimensions = true
					}
					aio.PrintMessage(fmt.Sprintf("Unsupported key %q in %q command, skipping", lit2, lit))
				}
				if err != nil {
					stopdata = true
				}
			}
		case FORMAT:
			// Format of the data bock: datatype, missing, gap
			stopformat := false
			for !stopformat {
				tok2, lit2 := p.scanIgnoreWhitespace()

				switch tok2 {
				case ENDOFCOMMAND:
					stopformat = true
				case DATATYPE:
					tok3, lit3 := p.scanIgnoreWhitespace()
					if tok3 != EQUAL {
						err = fmt.Errorf("Expecting '=' after DATATYPE, got %q", lit3)
						stopformat = true
					} else {
						tok4, lit4 := p.scanIgnoreWhitespace()
						if tok4 == IDENT {
							datatype = lit4
						} else {
							err = fmt.Errorf("Expecting identifier after 'DATATYPE=', got %q", lit4)
							stopformat = true
						}
					}
				case MISSING:
					tok3, lit3 := p.scanIgnoreWhitespace()
					if tok3 != EQUAL {
						err = fmt.Errorf("Expecting '=' after MISSING, got %q", lit3)
						stopformat = true
					} else {
						tok4, lit4 := p.scanIgnoreWhitespace()
						if tok4 != IDENT {
							err = fmt.Errorf("Expecting Integer value after 'MISSING=', got %q", lit4)
							stopformat = true
						} else {
							if len(lit4) != 1 {
								err = fmt.Errorf("Expecting a single character after MISSING=', got %q", lit4)
								stopformat = true
							} else {
								missing = []rune(lit4)[0]
							}
						}
					}
				case GAP:
					tok3, lit3 := p.scanIgnoreWhitespace()
					if tok3 != EQUAL {
						err = fmt.Errorf("Expecting '=' after GAP, got %q", lit3)
						stopformat = true
					} else {
						tok4, lit4 := p.scanIgnoreWhitespace()
						if tok4 != IDENT {
							err = fmt.Errorf("Expecting an identifier after 'GAP=', got %q", lit4)
							stopformat = true
						} else {
							if len(lit4) != 1 {
								err = fmt.Errorf("Expecting a single character after GAP=', got %q", lit4)
								stopformat = true
							} else {
								gap = []rune(lit4)[0]
							}
						}
					}
				case MATCHCHAR:
					tok3, lit3 := p.scanIgnoreWhitespace()
					if tok3 != EQUAL {
						err = fmt.Errorf("Expecting '=' after MATCHCHAR, got %q", lit3)
						stopformat = true
					} else {
						tok4, lit4 := p.scanIgnoreWhitespace()
						if tok4 != IDENT {
							err = fmt.Errorf("Expecting character value after 'MATCHCHAR=', got %q", lit4)
							stopformat = true
						} else {
							if len(lit4) != 1 {
								err = fmt.Errorf("Expecting a single character after MATCHCHAR=', got %q", lit4)
								stopformat = true
							} else {
								matchchar = []rune(lit4)[0]
							}
						}
					}
				default:
					if err = p.parseUnsupportedKey(lit2); err != nil {
						stopformat = true
					}
					aio.PrintMessage(fmt.Sprintf("Unsupported key %q in %q command, skipping", lit2, lit))
				}
				if err != nil {
					stopdata = true
				}
			}
		case MATRIX:
			// Character matrix (Alignmemnt)
			// So far: Does not handle interleave case...
			stopmatrix := false
			for !stopmatrix {
				tok2, lit2 := p.scanIgnoreWhitespace()
				switch tok2 {
				case OPENBRACK:
					if tok2, lit2, err = p.consumeComment(tok2, lit2); err != nil {
						stopmatrix = true
					}
				case IDENT, NUMERIC:
					// We remove whitespaces in sequences if any
					// and take into account possibly interleaved
					// sequences
					stopseq := false
					name := lit2
					sequence := ""
					for !stopseq {
						tok3, lit3 := p.scanIgnoreWhitespace()
						switch tok3 {
						case IDENT:
							sequence = sequence + lit3
						case ENDOFLINE:
							stopseq = true
						default:
							err = fmt.Errorf("Expecting sequence after sequence identifier (%q) in Matrix block, got %q", lit2, lit3)
							stopseq = true
						}
					}
					if err != nil {
						stopmatrix = true
					} else {
						addseq(sequences, &names, sequence, name)
					}
				case ENDOFLINE:
					break
				case ENDOFCOMMAND:
					stopmatrix = true
				default:
					err = fmt.Errorf(fmt.Sprintf("Expecting sequence identifier in Matrix block, got %q", lit2))
					stopmatrix = true
				}
			}
			if err != nil {
				stopdata = true
			}
		case OPENBRACK:
			if tok, lit, err = p.consumeComment(tok, lit); err != nil {
				stopdata = true
			}
		default:
			err = p.parseUnsupportedCommand()
			aio.PrintMessage(fmt.Sprintf("Unsupported command %q in block DATA, skipping", lit))
			if err != nil {
				stopdata = true
			}
		}
	}
	return
}

// Just skip the current command
func (p *Parser) parseUnsupportedCommand() (err error) {
	// Unsupported data command
	stopunsupported := false
	for !stopunsupported {
		tok, lit := p.scanIgnoreWhitespace()
		switch tok {
		case ILLEGAL:
			err = fmt.Errorf("found illegal token %q", lit)
			stopunsupported = true
		case EOF:
			err = fmt.Errorf("End of file within a command (no;)")
			stopunsupported = true
		case ENDOFCOMMAND:
			stopunsupported = true
		}
	}
	return
}

// Just skip the current key
func (p *Parser) parseUnsupportedKey(key string) (err error) {
	// Unsupported token
	tok, lit := p.scanIgnoreWhitespace()
	if tok != EQUAL {
		err = fmt.Errorf("Expecting '=' after %s, got %q", key, lit)
	} else {
		tok2, lit2 := p.scanIgnoreWhitespace()
		if tok2 != IDENT && tok2 != NUMERIC {
			err = fmt.Errorf("Expecting an identifier after '%s=', got %q", key, lit2)
		}
	}
	return
}

// Just skip the current block
func (p *Parser) parseUnsupportedBlock() error {
	var err error
	stopunsupported := false
	for !stopunsupported {
		tok, lit := p.scanIgnoreWhitespace()
		switch tok {
		case ILLEGAL:
			err = fmt.Errorf("found illegal token %q", lit)
			stopunsupported = true
		case EOF:
			err = fmt.Errorf("End of file within a block (no END;)")
			stopunsupported = true
		case END:
			tok2, _ := p.scanIgnoreWhitespace()
			if tok2 != ENDOFCOMMAND {
				err = fmt.Errorf("End token without ;")
			}
			stopunsupported = true
		}
	}
	return err
}

// Append a sequence to the hashmap map[name]seq.
// If the name does not exists in the map, adds the name to names and the sequence to the map (to keep order of the input nexus matrix
// Otherwise, append the sequence to the already insert sequence in the map
func addseq(sequences map[string]string, names *[]string, sequence string, name string) {
	seq, ok := sequences[name]
	if !ok {
		*names = append(*names, name)
		sequences[name] = sequence
	} else {
		sequences[name] = seq + sequence
	}
}

// Consumes comment inside brakets [comment] if the given current token is a [.
// At the end returns the matching ] token and lit.
// If the given token is not a [, then returns the input token and lit
func (p *Parser) consumeComment(curtoken Token, curlit string) (outtoken Token, outlit string, err error) {
	outtoken = curtoken
	curlit = curlit
	if curtoken == OPENBRACK {
		for outtoken != CLOSEBRACK {
			outtoken, outlit = p.scanIgnoreWhitespace()
			if outtoken == EOF || outtoken == ILLEGAL {
				err = fmt.Errorf("Unmatched bracket")
			}
		}
	}
	return
}
