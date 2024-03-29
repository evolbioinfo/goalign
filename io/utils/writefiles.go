package utils

import (
	"bufio"
	"compress/gzip"
	"io"
	"os"
	"strings"

	"github.com/ulikunitz/xz"
)

type StringWriterCloser interface {
	io.Writer
	io.Closer
	io.StringWriter
}

type gzstringwritercloser struct {
	f   *os.File
	gw  *gzip.Writer
	buf *bufio.Writer
}

type xzstringwritercloser struct {
	f   *os.File
	xw  *xz.Writer
	buf *bufio.Writer
}

func (gswc *gzstringwritercloser) Close() (err error) {
	if err = gswc.buf.Flush(); err != nil {
		return
	}
	if err = gswc.gw.Close(); err != nil {
		return
	}
	return gswc.f.Close()
}

func (gswc *gzstringwritercloser) Write(p []byte) (nn int, err error) {
	return gswc.buf.Write(p)
}

func (gswc *gzstringwritercloser) WriteString(s string) (nn int, err error) {
	return gswc.buf.WriteString(s)
}

func (xswc *xzstringwritercloser) Close() (err error) {
	if err = xswc.buf.Flush(); err != nil {
		return
	}
	if err = xswc.xw.Close(); err != nil {
		return
	}
	return xswc.f.Close()
}

func (xswc *xzstringwritercloser) Write(p []byte) (nn int, err error) {
	return xswc.buf.Write(p)
}

func (xswc *xzstringwritercloser) WriteString(s string) (nn int, err error) {
	return xswc.buf.WriteString(s)
}

func OpenWriteFile(file string) (f StringWriterCloser, err error) {
	if file == "stdout" || file == "-" {
		f = os.Stdout
	} else if file == "none" {
		f, err = os.OpenFile(os.DevNull, os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0666)
	} else if strings.HasSuffix(file, ".gz") {
		var fi *os.File
		if fi, err = os.Create(file); err != nil {
			return
		}
		gw := gzip.NewWriter(fi)
		buf := bufio.NewWriter(gw)
		f = &gzstringwritercloser{f: fi, gw: gw, buf: buf}
	} else if strings.HasSuffix(file, ".xz") {
		var fi *os.File
		if fi, err = os.Create(file); err != nil {
			return
		}
		xw, _ := xz.NewWriter(fi)
		buf := bufio.NewWriter(xw)
		f = &xzstringwritercloser{f: fi, xw: xw, buf: buf}
	} else {
		f, err = os.Create(file)
	}
	return
}

func CloseWriteFile(f io.Closer, filename string) {
	if filename != "-" && filename != "stdout" && filename != "none" {
		f.Close()
	}
}
