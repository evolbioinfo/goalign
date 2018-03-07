/*
Package intended to draw alignmentd on different devices :
 - Terminal,
 - Html file
 - ...
*/
package draw

import (
	"github.com/fredericlemoine/goalign/align"
)

/*
Generic struct that represents tree layout:
 * normal layout...
*/
type AlignLayout interface {
	DrawAlign(a align.Alignment) error
}
