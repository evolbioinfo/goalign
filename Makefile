GO_EXECUTABLE := go
VERSION := $(shell git describe --abbrev=10 --dirty --always --tags)
VERSION_PACKAGE := github.com/evolbioinfo/goalign/version.Version
NAME := goalign
PACKAGE:=github.com/evolbioinfo/goalign

all: dep build test 

dep:
	${GO_EXECUTABLE} get .

build:
	${GO_EXECUTABLE} build -o ${NAME} -ldflags "-X ${VERSION_PACKAGE}=${VERSION}" ${PACKAGE}

install: testcommands 
	rm -f ${GOPATH}/bin/${NAME}
	${GO_EXECUTABLE} install -ldflags "-X ${VERSION_PACKAGE}=${VERSION}" ${PACKAGE}

test: dep
	${GO_EXECUTABLE} test ${PACKAGE}/...

testcommands:
	bash ./test.sh

.PHONY: deploy deploydir deploywinamd deploywin386 deploylinuxamd deploylinux386 deploydarwinamd

deploy: deploywinamd deploywin386 deploylinuxamd deploylinux386 deploydarwinamd deploydarwinarm
	tar -czvf deploy/${VERSION}.tar.gz --directory="deploy" ${VERSION}

deploydir:
	mkdir -p deploy/${VERSION}

deploywinamd: dep deploydir
	env GOOS=windows GOARCH=amd64 ${GO_EXECUTABLE} build -o deploy/${VERSION}/${NAME}_${VERSION}_amd64.exe -ldflags "-X ${VERSION_PACKAGE}=${VERSION}" ${PACKAGE}

deploylinuxamd: dep deploydir
	env GOOS=linux GOARCH=amd64 ${GO_EXECUTABLE} build -o deploy/${VERSION}/${NAME}_${VERSION}_amd64_linux -ldflags "-X ${VERSION_PACKAGE}=${VERSION}" ${PACKAGE}

deploydarwinamd: dep deploydir
	env GOOS=darwin GOARCH=amd64 ${GO_EXECUTABLE} build -o deploy/${VERSION}/${NAME}_${VERSION}_amd64_darwin -ldflags "-X ${VERSION_PACKAGE}=${VERSION}" ${PACKAGE}

deploydarwinarm: dep deploydir
	env GOOS=darwin GOARCH=arm64 ${GO_EXECUTABLE} build -o deploy/${VERSION}/${NAME}_${VERSION}_arm64_darwin -ldflags "-X ${VERSION_PACKAGE}=${VERSION}" ${PACKAGE}
