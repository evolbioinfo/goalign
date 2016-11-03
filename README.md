# GoAlign
GoAlign is a set of command line tools to manipulate phylogenetic trees. It is implemented in [Go](https://golang.org/) language.

## Installation
### Binaries
You can download already compiled binaries for the latest release in the [release](https://github.com/fredericlemoine/goalign/releases) section.
Binaries are available for MacOS, Linux, and Windows (32 and 64 bits).

Once downloaded, you can just run the executable without any other downloads.

### From sources
In order to compile GoAlign, you must first [download](https://golang.org/dl/) and [install](https://golang.org/doc/install) Go on your system.

Then you just have to type :
```
go get github.com/fredericlemoine/goalign/
```
This will download GoAlign sources from github, and all its dependencies.

You can then build it with:
```
cd $GOPATH/src/github.com/fredericlemoine/goalign/
make
```
The `goalign` executable should be located in the `$GOPATH/bin` folder.

## Usage
