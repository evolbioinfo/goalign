package countprofile

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io/utils"
)

// FromFile Parses a "profile" file constisting of a number of occurences of each character
// per site, tab separated, in the form:
// site - A B C D G H K M N R...
// 0 1 2 3 4 0 ...
func FromFile(file string) (p *align.CountProfile, err error) {
	var f *os.File
	var r *bufio.Reader
	var gr *gzip.Reader
	var l string
	var i int
	var field string

	p = align.NewCountProfile()

	if file == "stdin" || file == "-" {
		f = os.Stdin
	} else {
		if f, err = os.Open(file); err != nil {
			return
		}
	}

	if strings.HasSuffix(file, ".gz") {
		if gr, err = gzip.NewReader(f); err != nil {
			return
		}
		r = bufio.NewReader(gr)
	} else {
		r = bufio.NewReader(f)
	}

	// We parse the header
	if l, err = utils.Readln(r); err != nil {
		return
	}
	headslice := strings.Split(l, "\t")
	header := make([]rune, 0, len(headslice)-1)
	for i, field = range headslice {
		if i > 0 {
			r := []rune(field)
			if len(r) != 1 {
				err = fmt.Errorf("Character name Should be One character")
				return
			}
			header = append(header, r[0])
		}
	}
	p.SetHeader(header)

	// Then the counts
	var count int = 0
	l, err = utils.Readln(r)
	for err == nil {
		for i, field = range strings.Split(l, "\t") {
			if i > 0 {
				if count, err = strconv.Atoi(field); err != nil {
					return
				}
				p.AppendCount(i-1, count)
			}
		}
		l, err = utils.Readln(r)
	}
	if err.Error() == "EOF" {
		err = nil
	}
	return
}
