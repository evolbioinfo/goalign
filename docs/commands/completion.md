# Goalign: toolkit and api for alignment manipulation

## Commands

### completion
This command generates auto-completion commands for bash or zsh. 

#### Usage

```
Usage:
  goalign completion SHELL
```

#### Bash
* Install bash-completion:
```
# MacOS brew
brew install bash-completion
# MacOS port (do not forget to change
# the path to bash command in terminal
# preferences to /opt/local/bin/bash -l)
sudo port install bash-completion
# Linux
yum install bash-completion -y
apt-get install bash-completion
```

* Activate goalign bash completion
```
# Once
source <(goalign completion bash)
# Permanently
mkdir ~/.goalign
goalign completion bash > ~/.goalign/completion.bash.inc
printf "
# goalign shell completion
source '$HOME/.goalign/completion.bash.inc'
" >> $HOME/.bashrc
```

#### Zsh (not tested)

```
# Once
source <(goalign completion zsh)
# Permanently
goalign completion zsh > "${fpath[1]}/_goalign"
```
