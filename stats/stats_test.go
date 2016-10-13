package stats

import (
	"fmt"
	"os"
	"testing"
)

func TestDirichlet(t *testing.T) {

	size := 1000
	alpha := make([]float64, size)
	for i := 0; i < size; i++ {
		alpha[i] = 1.0
	}
	sample := Dirichlet(alpha...)

	if len(sample) != size {
		t.Error("Size of sample is different from alpha slice")
	}

	for _, a := range sample {
		if a < 0 {
			t.Error("Dirichlet Sample should be positive")
		}
		fmt.Fprintf(os.Stdout, "\t%f", a)
	}
}
