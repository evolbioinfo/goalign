package align

import (
	"fmt"
)

// CountProfile represents a simple view of an alignment
// and stores the number of occurences of each characters at
// each position of an alignment
type CountProfile struct {
	names  []int   // mapping between character and index. names uses utf8 int value to store indexes names[int(r)]=index
	header []uint8 // characters in the count table, similar to prev map
	counts [][]int // first dimension: len(head), second dimension; nb sites in the underlying alignment
}

// NewCountProfileFromAlignment initializes a new CountProfile using an input alignment
func NewCountProfileFromAlignment(al Alignment) (p *CountProfile) {
	p = NewCountProfile()

	p.header = make([]uint8, 0, 100)
	p.counts = make([][]int, 0, 100)

	al.IterateChar(func(name string, seq []uint8) bool {
		for i, r := range seq {
			idx := p.names[int(r)]
			if idx < 0 {
				idx = len(p.header)
				p.names[int(r)] = idx
				p.header = append(p.header, r)
				p.counts = append(p.counts, make([]int, al.Length()))
			}
			p.counts[idx][i]++
		}
		return false
	})
	return
}

// NewCountProfile initializes a new Profile with nil attributes
func NewCountProfile() (p *CountProfile) {
	p = &CountProfile{
		names:  make([]int, 130),
		header: nil,
		counts: nil,
	}

	// Initialize to "not exit" all characters
	for i := 0; i < 130; i++ {
		p.names[i] = -1
	}

	return
}

// CountsAt returns the counts for all sites, for the ith character
// (arbitrary order of character)
func (p *CountProfile) CountsAt(i int) (counts []int, err error) {
	if i > len(p.counts) || i < 0 {
		err = fmt.Errorf("No counts for character at index %d", i)
		return
	}
	counts = p.counts[i]
	return
}

// NameAt returns the name of ith character in the header
func (p *CountProfile) NameAt(i int) (name uint8, err error) {
	if i >= len(p.header) || i < 0 {
		err = fmt.Errorf("No character at index %d", i)
		return
	}
	name = p.header[i]
	return
}

// NameIndex returns the index of the given character in the header
// If the character does not exist, returns false
func (p *CountProfile) NameIndex(r uint8) (index int, ok bool) {
	index = p.names[int(r)]
	ok = (index >= 0)
	return
}

// Count returns the number of occurences of the character r at the position site
func (p *CountProfile) Count(r uint8, site int) (count int, err error) {
	var index int
	count = 0
	if index = p.names[int(r)]; index < 0 {
		err = fmt.Errorf("Character does not exist in the Profile %c", r)
		return
	}
	if site < 0 || site >= len(p.counts[index]) {
		err = fmt.Errorf("Site %d is outside the profile", site)
		return
	}
	count = p.counts[index][site]
	return
}

// CountAt returns the number of occurences of the ith character at the position site
func (p *CountProfile) CountAt(i, site int) (count int, err error) {
	if i >= len(p.header) || i < 0 {
		err = fmt.Errorf("No character at index %d", i)
		return
	}
	if site < 0 || site >= len(p.counts[i]) {
		err = fmt.Errorf("Site %d is outside the profile", site)
		return
	}
	count = p.counts[i][site]
	return
}

// NbCharacters returns the number of different characters in the profile
func (p *CountProfile) NbCharacters() (nb int) {
	return len(p.header)
}

// CheckLength returns true if the number of sites of the profile corresponds to the given length
// false otherwise.
func (p *CountProfile) CheckLength(length int) bool {
	for _, v := range p.counts {
		if len(v) != length {
			return false
		}
	}
	return true
}

// SetHeader sets the Header and initializes the count structure
func (p *CountProfile) SetHeader(header []uint8) {
	p.header = header
	p.counts = make([][]int, len(p.header))
	for i, r := range header {
		p.names[int(r)] = i
		p.counts[i] = make([]int, 0, 100)
	}
	return
}

// AppendCount appends a new site to the profile for the ith character, and
// associates count to it
func (p *CountProfile) AppendCount(i, count int) (err error) {
	if i >= len(p.header) || i < 0 {
		err = fmt.Errorf("No character at index %d", i)
		return
	}
	p.counts[i] = append(p.counts[i], count)
	return
}

func (p *CountProfile) Print() {
	fmt.Print("site")
	for _, n := range p.header {
		fmt.Printf("\t%c", n)
	}
	fmt.Println()
	i := 0
	length := 0
	for {
		fmt.Printf("%d", i)
		for _, c := range p.counts {
			length = len(c)
			fmt.Printf("\t%d", c[i])
		}
		fmt.Println()
		i++
		if i >= length {
			break
		}
	}
}
