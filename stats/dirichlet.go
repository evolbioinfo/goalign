package stats

import (
	"errors"
	"fmt"
	"github.com/evolbioinfo/goalign/io"
)

// Dirichlet returns a set of random numbers from dirichlet distribution
// See https://en.wikipedia.org/wiki/Dirichlet_distribution#Random_number_generation
// For further information
// Factor: multiplying factor to apply to all values : Normally it is 1
func Dirichlet(factor float64, alpha ...float64) []float64 {
	if len(alpha) <= 2 {
		io.ExitWithMessage(errors.New(fmt.Sprintf("Alpha parameter vector for Dirichlet sample contains less than 2 values")))
	}
	sample := make([]float64, len(alpha))
	sum := 0.0
	for i, a := range alpha {
		if a <= 0.0 {
			io.ExitWithMessage(errors.New(fmt.Sprintf("Invalid parameter alpha %.2f", a)))
		}
		sample[i] = gamma(a, 1)
		sum += sample[i]
	}
	for i, _ := range alpha {
		sample[i] = factor * sample[i] / sum
	}
	return sample
}
