package stats

import (
	"fmt"
	"math"
	"os"
	"testing"
)

func TestDirichlet(t *testing.T) {

	size := 1000
	alpha := make([]float64, size)
	for i := range size {
		alpha[i] = 1.0
	}
	sample, _ := Dirichlet(1.0, alpha...)

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

func TestDirichlet1(t *testing.T) {

	size := 1000
	factor := float64(size)
	sample, _ := Dirichlet1(factor, size)

	if len(sample) != size {
		t.Error("Size of sample is different from alpha slice")
	}

	sum := 0.0
	for _, a := range sample {
		sum += a
		if a < 0 {
			t.Error("Dirichlet Sample should be positive")
		}
		fmt.Fprintf(os.Stdout, "\t%f", a)
	}
	if math.Abs(sum-factor) > 0.00000000001 {
		t.Errorf("Dirichlet sum %f != %f", sum, factor)
	}
}
