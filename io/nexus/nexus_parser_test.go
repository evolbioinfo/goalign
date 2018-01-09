package nexus_test

import (
	"strings"
	"testing"

	"github.com/fredericlemoine/goalign/io/nexus"
)

// Ensure the parser can parse strings into Statement ASTs.
func TestParser_ParseTree(t *testing.T) {
	goodnexus := [...]string{
		`#NEXUS
BEGIN TAXA;
      TaxLabels fish frog snake mouse;
END;

BEGIN CHARACTERS;
      Dimensions NChar=40;
      Format DataType=DNA;
      Matrix
        fish   ACATA GAGGG TACCT CTAAA
        fish   ACATA GAGGG TACCT CTAAG

        frog   ACATA GAGGG TACCT CTAAC
        frog   CCATA GAGGG TACCT CTAAG

        snake  ACATA GAGGG TACCT CTAAG
        snake  GCATA GAGGG TACCT CTAAG

        mouse  ACATA GAGGG TACCT CTAAT
        mouse  TCATA GAGGG TACCT CTAAG
;
END;
EOF
`,
	}

	for i, innexus := range goodnexus {
		align, err := nexus.NewParser(strings.NewReader(innexus)).Parse()
		if err != nil {
			t.Errorf("Nexus parser error %d ERROR: %s\n", i, err.Error())
		} else {
			if align.NbSequences() != 4 {
				t.Errorf("Alignment should have 4 sequences but has %d\n", align.NbSequences())
			}
			if align.Length() != 40 {
				t.Errorf("Alignment should be 40 nt long, but is %d\n", align.Length())
			}
		}
	}
}
