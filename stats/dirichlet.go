package stats

import (
	"fmt"
	"math/rand"
	"slices"
)

// Dirichlet returns a set of random numbers from dirichlet distribution
// (See https://en.wikipedia.org/wiki/Dirichlet_distribution#Random_variate_generation)
// Factor: multiplying factor to apply to all values
func Dirichlet(factor float64, alpha ...float64) (sample []float64, err error) {
	if len(alpha) <= 2 {
		err = fmt.Errorf("alpha parameter vector for Dirichlet sample contains less than 2 values")
		return
	}
	sample = make([]float64, len(alpha))
	sum := 0.0
	for i, a := range alpha {
		if a <= 0.0 {
			err = fmt.Errorf("invalid parameter alpha %.2f", a)
			return
		}
		sample[i] = gamma(a, 1)
		sum += sample[i]
	}
	for i, _ := range alpha {
		sample[i] = factor * sample[i] / sum
	}
	return
}

// Dirichlet1 returns a set of random numbers from dirichlet distribution
// When each alpha equals 1 (see https://en.wikipedia.org/wiki/Dirichlet_distribution#When_each_alpha_is_1)
func Dirichlet1(factor float64, nvalues int) (sample []float64, err error) {
	if nvalues <= 2 {
		err = fmt.Errorf("nvalues should be > 2")
		return
	}
	sample = make([]float64, nvalues)
	intervals := make([]float64, nvalues+1)
	intervals[0] = 0.0
	intervals[1] = 1.0
	for i := 2; i < nvalues+1; i++ {
		intervals[i] = rand.Float64()
	}
	slices.Sort(intervals)
	for i := 1; i < nvalues+1; i++ {
		sample[i-1] = factor * (intervals[i] - intervals[i-1])
	}
	return
}
