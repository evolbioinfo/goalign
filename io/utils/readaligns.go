package utils

import (
	"bufio"
	"errors"
	"io"

	"github.com/fredericlemoine/goalign/align"
	"github.com/fredericlemoine/goalign/io/fasta"
	"github.com/fredericlemoine/goalign/io/phylip"
	tnexus "github.com/fredericlemoine/gotree/io/nexus"
)

// Parses the input buffer while automatically
// detecting the format between Newick, Phylip, and Nexus
//
// If several alignments are present in the onput file, only the first will be
// parsed.
//
// Returned format may be align.FORMAT_PHYLIP, align.FORMAT_FASTA, or align.FORMAT_NEXUS
//
// rootinpustrict: In the case of phylip detected format: should we consider it as strict or not?
//
// There is no new go routine here because only 1 alignment is parsed. No need to give a closer.
// If the reader comes from a file, the file must be closed in the calling function.
func ParseAlignmentAuto(r *bufio.Reader, rootinputstrict bool) (al align.Alignment, format int, err error) {
	var firstbyte byte
	var nex *tnexus.Nexus

	if firstbyte, err = r.ReadByte(); err != nil {
		return
	}

	if err = r.UnreadByte(); err != nil {
		return
	}

	// First test Fasta format
	if firstbyte == '>' {
		format = align.FORMAT_FASTA
		al, err = fasta.NewParser(r).Parse()
	} else if firstbyte == '#' {
		if nex, err = tnexus.NewParser(r).Parse(); err != nil {
			return
		}
		format = align.FORMAT_NEXUS
		// Then test nexus
		if !nex.HasAlignment {
			err = errors.New("Nexus format has no alignment")
		}
		al = nex.Alignment()
	} else {
		// Finally test Phylip
		format = align.FORMAT_PHYLIP
		al, err = phylip.NewParser(r, rootinputstrict).Parse()
	}

	return
}

// Parses the input buffer while automatically
// detecting the format between Newick, Phylip, and Nexus
//
// If several alignments are present in the input file, they are queued in the channel
//
// rootinpustrict: In the case of phylip detected format: should we consider it as strict or not?
//
// If there is something to close ( f!=nil) after the parsing (like input file, etc.), f will be closed
// after parsing is finished (even in the go routine in the case of several input alignments).
// If the alignment comes from a file for exemple, the file will be closed by this function, so no need to
// do it in the calling function
func ParseMultiAlignmentsAuto(f io.Closer, r *bufio.Reader, rootinputstrict bool) (alchan align.AlignChannel, format int, err error) {
	var al align.Alignment
	alchan.Achan = make(chan align.Alignment, 15)
	var firstbyte byte
	var nex *tnexus.Nexus

	if firstbyte, err = r.ReadByte(); err != nil {
		return
	}

	if err = r.UnreadByte(); err != nil {
		return
	}
	// First test Fasta format
	if firstbyte == '>' {
		if al, err = fasta.NewParser(r).Parse(); err != nil {
			return
		}
		format = align.FORMAT_FASTA
		alchan.Achan <- al
		if f != nil {
			f.Close()
		}
		close(alchan.Achan)
	} else if firstbyte == '#' {
		if nex, err = tnexus.NewParser(r).Parse(); err != nil {
			return
		}
		format = align.FORMAT_NEXUS
		// Then test nexus
		if !nex.HasAlignment {
			return
		}
		alchan.Achan <- nex.Alignment()
		if f != nil {
			f.Close()
		}
		close(alchan.Achan)
	} else {
		format = align.FORMAT_PHYLIP
		// Finally test Phylip
		go func() {
			if err := phylip.NewParser(r, rootinputstrict).ParseMultiple(alchan.Achan); err != nil {
				alchan.Err = err
			}
			if f != nil {
				f.Close()
			}
		}()
	}
	return
}
