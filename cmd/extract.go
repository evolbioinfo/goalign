package cmd

import (
	"bufio"
	"compress/gzip"
	"errors"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/evolbioinfo/goalign/align"
	"github.com/evolbioinfo/goalign/io"

	"github.com/evolbioinfo/goalign/io/utils"
	"github.com/spf13/cobra"
)

type extractSubSequence struct {
	starts []int
	ends   []int
	name   string
	strand bool
}

var extractrefseq string
var extractcoordfile string
var extractoutput string
var extracttranslate int
var extractgff bool
var extractsuffix string
var extractprefix string

// extractCmd represents the extract command
var extractCmd = &cobra.Command{
	Use:   "extract",
	Short: "Extracts several sub-sequences from an input alignment",
	Long: `This command extracts several sub-alignments from an input alignment. 
	
	It is similar to subseq.md, with two main differences:

	1. As input, it takes an annotation file defining the coordinates of sub-alignments 
	   to extract, and can then extract several sub-alignments in one command;
	2. Each sub-alignment may be defined as several "blocks" (~bed format blocks), 
	   potentially overlapping (even in any order).
	
	This command takes an alignment and extracts several sub-alignments from it. 
	Subs-alignments are defined in an input tab separated file with the following mandatory 
	columns:
	
	1. Start coordinates: 0-based inclusive. If the sub-alignment is defined by several blocks, 
	   several start coordinates may be given, and separated by comas;
	2. End coordinates: 1-based (or 0-based exclusive). If the sub-alignment is defined by several 
	   blocks, several end coordinates may be given (same number as start coordinates), and 
	coma separated;
	3. Name of the subsequence
	
	Example of an annotation file:
	
	0	10	orf1
	10	100	orf2
	100,105	106,110	orf3
	
	The 3rd line defines a sub-alignment containing positions [100-106[+[105-110[ (or [100-105]+[105-109]).
	
	If start (or end) coordinates are outside the alignment, or are not compatible (start>=end)
	then it exits with an error.
	
	If a sub-alignment is defined by several blocks, they are allowed to overlap or be in any order.
	
	Output file names will be defined by the names of the subsequences, with .fa or .phy 
	extension depending on the input file format.
	
	If --ref-seq is given, then the coordinates are defined wrt the given reference sequence 
	(gaps are not taken into acount although they are still present in the output sub-alignment).
	
	If --translate is >=0 and the input alignment is nucleotidic, extracted subsequences are translated 
	into amino acids.
	- If --translate < 0 : No translation
	- If --translate 0: Standard genetic code
	- If --translate 1: Vertebrate mitochondrial genetic code
	- If --translate 2: Invertebrate mitochondrial genetic code
	
	If --ref-seq is given, and --translate>=0 is given, be careful about the behavior! It will extract the 
	nucleotide sequences corresponding to start/stop coordinates on the reference sequence, and will translate
	the resulting alignment with the following process: The extracted alignment will be translated codon by
	codon using the given reference sequence as guide, by iterating over the reference non gap nucleotides 3 by 3. 
	At each iteration, the current reference codon may have gaps between nucleotides, and the translation of the
	current codon will be done as following:
		* ex 1:
			Ref: AC--GTACGT
			Seq: ACTTGTACGT
			In that case, the first ref codon is [0,1,4], corresponding to sequence ACTTG in seq
			ACTTG % 3 != 0 ==> Frameshift? => Replaced by T in ref and X in the compared sequence.
		* ex 2:
			Ref: AC---GTACGT
			Seq: ACTTTGTACGT
			ref codon: [0,1,5]
			seq      : ACTTTG (%3==0): Insertion - OK => Replaced by "T-" in ref and "TL" in seq
		* ex 3:
			Ref: ACGTACGT
			Seq: A--TACGT
			ref codon: [0,1,2]
			seq      : A--: Deletion: not ok : Frameshift? => Replaced by "T" in ref and "X" in comp
		* ex 4:
			Ref: AC----GTACGT
			Seq: ACTT-TGTACGT
			ref codon: [0,1,6]
			seq      : ACTTTG (%3==0): Insertion - OK => Replaced by "T-" in ref and "TT" in seq
		* ex 5:
			Ref: AC----GTACGT
			Seq: ACT--TGTACGT
			ref codon: [0,1,6]
			seq      : ACTTTG : Insertion not OK : Frameshift? => Replaced by "T-" in ref and "XX" in seq
	This allows to easily translate a multiple sequence alignment containing partial sequences, but the 
	interpretation should be careful: the translation of some sequences may not be representative of the 
	translation of the unaligned sequences.


	Basic example:
	goalign extract -i alignment.fasta -f annotations.txt
	
	If the input file contains several alignments, only the first one is considered.
`,
	RunE: func(cmd *cobra.Command, args []string) (err error) {
		var aligns *align.AlignChannel
		var f utils.StringWriterCloser
		var subcoords []extractSubSequence
		var subalign, subaligntmp align.Alignment

		if extractcoordfile == "none" {
			err = fmt.Errorf("subsequence coordinate file should be specified")
			return
		}

		if aligns, err = readalign(infile); err != nil {
			io.LogError(err)
			return
		}

		if extractgff {
			if subcoords, err = parseGFFFile(extractcoordfile); err != nil {
				io.LogError(err)
				return
			}
		} else {
			if subcoords, err = parseCoordinateFile(extractcoordfile); err != nil {
				io.LogError(err)
				return
			}
		}

		refseq := cmd.Flags().Changed("ref-seq")

		al := <-aligns.Achan
		if aligns.Err != nil {
			err = aligns.Err
			io.LogError(err)
			return
		}

		for _, subseq := range subcoords {
			subalign = nil
			for i, s := range subseq.starts {
				e := subseq.ends[i]
				l := e - s
				if s < 0 || e > al.Length() {
					err = fmt.Errorf("coordinates are outside alignment: [%d,%d[", s, e)
					io.LogError(err)
					return
				}
				if s >= e {
					err = fmt.Errorf("block length should be >0 : [%d,%d[", s, e)
					io.LogError(err)
					return
				}

				if refseq {
					if s, l, err = al.RefCoordinates(extractrefseq, s, l); err != nil {
						io.LogError(err)
						return
					}
				}
				if subaligntmp, err = al.SubAlign(s, l); err != nil {
					io.LogError(err)
					return
				}

				if subalign == nil {
					subalign = subaligntmp
				} else {
					subalign.Concat(subaligntmp)
				}
			}

			if !subseq.strand {
				if err = subalign.ReverseComplement(); err != nil {
					io.LogError(err)
					return
				}
			}

			if al.Alphabet() == align.NUCLEOTIDS && extracttranslate >= 0 {
				if refseq {
					if err = subalign.TranslateByReference(0, extracttranslate, extractrefseq); err != nil {
						io.LogError(err)
						return
					}
				} else {
					if err = subalign.Translate(0, extracttranslate); err != nil {
						io.LogError(err)
						return
					}
				}
			}
			if f, err = utils.OpenWriteFile(fmt.Sprintf("%s%c%s%s%s%s", extractoutput, os.PathSeparator, extractprefix, subseq.name, extractsuffix, alignExtension())); err != nil {
				io.LogError(err)
				return
			}
			writeAlign(subalign, f)
			f.Close()
		}

		return
	},
}

func init() {
	RootCmd.AddCommand(extractCmd)
	extractCmd.PersistentFlags().StringVar(&extractrefseq, "ref-seq", "none", "Reference sequence on which coordinates are given")
	extractCmd.PersistentFlags().IntVar(&extracttranslate, "translate", -1, "Wether the extracted sequence will be translated (only if input alignment is nucleotide). <0: No translation, 0: Std code, 1: Vertebrate mito, 2: Invertebrate mito")
	extractCmd.PersistentFlags().StringVarP(&extractoutput, "output", "o", ".", "Output folder")
	extractCmd.PersistentFlags().StringVar(&extractcoordfile, "coordinates", "none", "File with all coordinates of the sequences to extract")
	extractCmd.PersistentFlags().BoolVar(&extractgff, "gff", false, "Wether the coordinate file specified with --coordinates is in gff format")
	extractCmd.PersistentFlags().StringVar(&extractprefix, "prefix", "", "The prefix of the generated files (before gene name)")
	extractCmd.PersistentFlags().StringVar(&extractsuffix, "suffix", "", "The suffix of the generated files (before extension)")

}

func parseCoordinateFile(file string) (coords []extractSubSequence, err error) {
	var f *os.File
	var r *bufio.Reader
	var gr *gzip.Reader
	var si int
	var strand bool

	coords = make([]extractSubSequence, 0)

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
	l, e := utils.Readln(r)
	for e == nil {
		cols := strings.Split(l, "\t")
		if cols == nil || (len(cols) != 3 && len(cols) != 4) {
			err = errors.New("bad format from coordinate file: There should be 3 or 4 columns")
			return
		}

		strand = true
		if len(cols) == 4 {
			strand = (cols[3] == "+")
		}
		subseq := extractSubSequence{
			starts: make([]int, 0),
			ends:   make([]int, 0),
			name:   cols[2],
			strand: strand,
		}

		startstr := strings.Split(cols[0], ",")
		endstr := strings.Split(cols[1], ",")

		if len(startstr) == 0 || len(endstr) == 0 || len(startstr) != len(endstr) {
			err = errors.New("bad format from coordinate file: start en end coordinates should have at least 1 coordinate and have the same length")
			return
		}

		for i, s := range startstr {
			if si, err = strconv.Atoi(s); err != nil {
				return
			}
			subseq.starts = append(subseq.starts, si)
			if si, err = strconv.Atoi(endstr[i]); err != nil {
				return
			}
			subseq.ends = append(subseq.ends, si)
		}
		coords = append(coords, subseq)
		l, e = utils.Readln(r)
	}

	return
}

func parseGFFFile(file string) (coords []extractSubSequence, err error) {
	var f *os.File
	var r *bufio.Reader
	var gr *gzip.Reader
	var ok bool
	var strand bool
	var start, end int
	var curgene *extractSubSequence
	// Name of the gene to the gene coords
	var coordsMap map[string]*extractSubSequence
	// Mapping between gene ids and gene names
	var geneID2Name map[string]string

	coordsMap = make(map[string]*extractSubSequence)
	geneID2Name = make(map[string]string)

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
	l, e := utils.Readln(r)

	for e == nil {
		cols := strings.Split(l, "\t")
		if cols == nil || len(cols) < 9 {
			err = errors.New("bad format from gff file: There should be 9 columns")
			return
		}

		infos := strings.Split(cols[8], ";")
		strand = cols[6] == "+"
		if start, err = strconv.Atoi(cols[3]); err != nil {
			err = fmt.Errorf("cannot convert start to int")
			return
		}
		if end, err = strconv.Atoi(cols[4]); err != nil {
			err = fmt.Errorf("cannot convert end to int")
			return
		}

		// We extract gene id and name from gene line
		if cols[2] == "gene" {
			gname := ""
			gid := ""
			for _, info := range infos {
				cols := strings.Split(info, "=")
				if len(cols) != 2 {
					err = fmt.Errorf("info field is malformed, a field is not in the form key=value")
					return
				}
				if cols[0] == "ID" {
					gid = cols[1]
				}
				if cols[0] == "Name" {
					gname = cols[1]
				}
			}
			geneID2Name[gid] = gname
		}

		// Now we want only want CDS lines
		if cols[2] == "CDS" {

			parentGeneID := ""
			parentGeneName := ""
			for _, info := range infos {
				cols := strings.Split(info, "=")
				if len(cols) != 2 {
					err = fmt.Errorf("info field is malformed, a field is not in the form key=value")
					return
				}
				if cols[0] == "Parent" {
					parentGeneID = cols[1]
				}
			}
			if parentGeneID == "" {
				err = fmt.Errorf("cds has no parent gene")
				return
			}
			if parentGeneName, ok = geneID2Name[parentGeneID]; !ok {
				err = fmt.Errorf("cds parent gene id has no name")
				return
			}

			if curgene, ok = coordsMap[parentGeneName]; !ok {
				curgene = &extractSubSequence{
					starts: make([]int, 0),
					ends:   make([]int, 0),
					name:   parentGeneName,
					strand: strand,
				}
				coordsMap[parentGeneName] = curgene
			}
			curgene.starts = append(curgene.starts, start-1)
			curgene.ends = append(curgene.ends, end)
		}
		l, e = utils.Readln(r)
	}

	for _, coord := range coordsMap {
		coords = append(coords, *coord)
	}

	return
}
