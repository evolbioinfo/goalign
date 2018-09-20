package cmd

import (
	"os"

	"github.com/spf13/cobra"
)

// itolCmd represents the itol command
var autocompleteCmd = &cobra.Command{
	Use:   "completion SHELL",
	Args:  cobra.ExactArgs(1),
	Short: "Generates auto-completion commands for bash or zsh",
	Long: `Generates auto-completion commands for bash or zsh. 

Examples (Largely inspired from kubectl command):
  # bash completion on macOS using homebrew
  ## If running Bash 3.2 included with macOS
  brew install bash-completion
  ## or, if running Bash 4.1+
  brew install bash-completion@2
  # Then add auto completion commands
  goalign completion bash > $(brew --prefix)/etc/bash_completion.d/goalign

  # bash completion on macOS using macports
  ## If running Bash 3.2 included with macOS
  Change path to bash in Terminal preferences/path to shell:
  /opt/local/bin/bash -l
  ## Then in stall bash-completion
  sudo port install bash-completion
  # Finally add auto completion commands
  goalign completion bash > $(brew --prefix)/etc/bash_completion.d/goalign

  # bash completion on Linux
  ## Load the goalign completion code for bash into the current shell
  source <(goalign completion bash)
  ## Write bash completion code to a file and source if from .bash_profile
  mkdir ~/.goalign
  goalign completion bash > ~/.goalign/completion.bash.inc
  printf "
  # goalign shell completion
  source '$HOME/.goalign/completion.bash.inc'
  " >> $HOME/.bashrc
  source $HOME/.bashrc

  # Load the goalign completion code for zsh[1] into the current shell
  source <(goalign completion zsh)
  # Set the goalign completion code for zsh[1] to autoload on startup
  goalign completion zsh > "${fpath[1]}/_goalign"

`,
	Run: func(cmd *cobra.Command, args []string) {
		shell := args[0]
		if shell == "bash" {
			RootCmd.GenBashCompletion(os.Stdout)
		} else if shell == "zsh" {
			RootCmd.GenZshCompletion(os.Stdout)
		}
	},
}

func init() {
	RootCmd.AddCommand(autocompleteCmd)
}
