package stats

// Inspired from https://github.com/leesper/go_rng/blob/master/gamma.go

import (
	"errors"
	"fmt"
	"math"
	"math/rand"

	"github.com/evolbioinfo/goalign/io"
)

// Gamma returns a random number of gamma distribution (alpha > 0.0 and beta > 0.0)
func Gamma(alpha, beta float64) float64 {
	if !(alpha > 0.0) || !(beta > 0.0) {
		io.ExitWithMessage(errors.New(fmt.Sprintf("Invalid parameter alpha %.2f beta %.2f", alpha, beta)))
	}
	return gamma(alpha, beta)
}

// inspired by random.py
func gamma(alpha, beta float64) float64 {
	var MAGIC_CONST float64 = 4 * math.Exp(-0.5) / math.Sqrt(2.0)
	if alpha > 1.0 {
		// Use R.C.H Cheng "The generation of Gamma variables with
		// non-integral shape parameters", Applied Statistics, (1977), 26, No. 1, p71-74

		ainv := math.Sqrt(2.0*alpha - 1.0)
		bbb := alpha - math.Log(4.0)
		ccc := alpha + ainv

		for {
			u1 := rand.Float64()
			if !(1e-7 < u1 && u1 < .9999999) {
				continue
			}
			u2 := 1.0 - rand.Float64()
			v := math.Log(u1/(1.0-u1)) / ainv
			x := alpha * math.Exp(v)
			z := u1 * u1 * u2
			r := bbb + ccc*v - x
			if r+MAGIC_CONST-4.5*z >= 0.0 || r >= math.Log(z) {
				return x * beta
			}
		}
	} else if alpha == 1.0 {
		u := rand.Float64()
		for u <= 1e-7 {
			u = rand.Float64()
		}
		return -math.Log(u) * beta
	} else { // alpha between 0.0 and 1.0 (exclusive)
		// Uses Algorithm of Statistical Computing - kennedy & Gentle
		var x float64
		for {
			u := rand.Float64()
			b := (math.E + alpha) / math.E
			p := b * u
			if p <= 1.0 {
				x = math.Pow(p, 1.0/alpha)
			} else {
				x = -math.Log((b - p) / alpha)
			}
			u1 := rand.Float64()
			if p > 1.0 {
				if u1 <= math.Pow(x, alpha-1.0) {
					break
				}
			} else if u1 <= math.Exp(-x) {
				break
			}
		}
		return x * beta
	}
}
