package cmd

import (
	"github.com/spf13/cobra"
)

// buildCmd represents the build command
var buildCmd = &cobra.Command{
	Use:   "build",
	Short: "Command to build bootstrap replicates",
	Long: `This command builds bootstrap replicates from an input alignment (fasta or phylip):

1. goalign build seqboot : Builds bootstrap alignments from an input alignment (nt or aa). Sequence order may be shuffled with option -S. Output alignments may be written in compressed files (--gz) and/or added in a tar archive (--tar).
2. goalign build distboot: Builds bootstrap distance matrices based on different models, from an input alignment (nt only). It builds n bootstrap alignments and computes a distance matrix for each replicate. All distance matrices are written in the output file. If the input alignment file contains several alignments, it will take the first one only. The following models for distance computation are available:
    - pdist
    - jc   : Juke-Cantor
    - k2p  : Kimura 2 Parameters
    - f81  : Felsenstein 81
    - f84  : Felsenstein 84
    - tn93 : Tamura and Nei 1993
`,
}

func init() {
	RootCmd.AddCommand(buildCmd)
}
