GO_EXECUTABLE := go
VERSION := $(shell git describe --abbrev=10 --dirty --always --tags)
DIST_DIRS := find * -type d -exec

all: build install

build:
	${GO_EXECUTABLE} build -o goalign -ldflags "-X github.com/fredericlemoine/goalign/cmd.Version=${VERSION}" github.com/fredericlemoine/goalign

install:
	${GO_EXECUTABLE} install -ldflags "-X github.com/fredericlemoine/goalign/cmd.Version=${VERSION}" github.com/fredericlemoine/goalign

test:
	${GO_EXECUTABLE} test github.com/fredericlemoine/goalign/tests/

deploy:
	mkdir -p deploy/${VERSION}
	env GOOS=windows GOARCH=amd64 ${GO_EXECUTABLE} build -o deploy/${VERSION}/goalign_amd64.exe -ldflags "-X github.com/fredericlemoine/goalign/cmd.Version=${VERSION}" github.com/fredericlemoine/goalign
	env GOOS=windows GOARCH=386 ${GO_EXECUTABLE} build -o deploy/${VERSION}/goalign_386.exe -ldflags "-X github.com/fredericlemoine/goalign/cmd.Version=${VERSION}" github.com/fredericlemoine/goalign
	env GOOS=darwin GOARCH=amd64 ${GO_EXECUTABLE} build -o deploy/${VERSION}/goalign_amd64_darwin -ldflags "-X github.com/fredericlemoine/goalign/cmd.Version=${VERSION}" github.com/fredericlemoine/goalign
	env GOOS=darwin GOARCH=386 ${GO_EXECUTABLE} build -o deploy/${VERSION}/goalign_386_darwin -ldflags "-X github.com/fredericlemoine/goalign/cmd.Version=${VERSION}" github.com/fredericlemoine/goalign
	env GOOS=linux GOARCH=amd64 ${GO_EXECUTABLE} build -o deploy/${VERSION}/goalign_amd64_linux -ldflags "-X github.com/fredericlemoine/goalign/cmd.Version=${VERSION}" github.com/fredericlemoine/goalign
	env GOOS=linux GOARCH=386 ${GO_EXECUTABLE} build -o deploy/${VERSION}/goalign_386_linux -ldflags "-X github.com/fredericlemoine/goalign/cmd.Version=${VERSION}" github.com/fredericlemoine/goalign
	tar -czvf deploy/${VERSION}.tar.gz deploy/${VERSION}
