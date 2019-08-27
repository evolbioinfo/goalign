package models

import (
	"math"
	"math/rand"

	"gonum.org/v1/gonum/stat/distuv"
)

const (
	DBL_MIN = 2.2250738585072014e-308
)

func IncompleteGamma(x, alpha, ln_gamma_alpha float64) float64 {
	/* returns the incomplete gamma ratio I(x,alpha) where x is the upper
	   	   limit of the integration and alpha is the shape parameter.
		      returns (-1) if in error
		      ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
		      (1) series expansion     if (alpha>x || x<=1)
		      (2) continued fraction   otherwise
		      RATNEST FORTRAN by
		      Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
		      19: 285-287 (AS32)
	              https://www.jstor.org/stable/2346339
	*/
	var i int
	var factor float64

	p := alpha
	g := ln_gamma_alpha
	accurate := 1.0e-8
	overflow := 1.0e30

	gin := 0.0
	rn := 0.0
	a := 0.0
	b := 0.0
	an := 0.0
	dif := 0.0
	term := 0.0
	pn := make([]float64, 6)

	if math.Abs(x) < DBL_MIN {
		return 0.0
	}

	if x < 0 || p <= 0 {
		return -1.0
	}

	factor = math.Exp(p*math.Log(x) - x - g)

	if x > 1 && x >= p {
		goto l30
	}
	/* (1) series expansion */
	gin = 1
	term = 1
	rn = p
l20:
	rn++
	term *= x / rn
	gin += term

	if term > accurate {
		goto l20
	}
	gin *= factor / p
	goto l50
l30:
	/* (2) continued fraction */
	a = 1 - p
	b = a + x + 1
	term = 0
	pn[0] = 1
	pn[1] = x
	pn[2] = x + 1
	pn[3] = x * b
	gin = pn[2] / pn[3]
l32:
	a++
	b += 2
	term++
	an = a * term
	for i = 0; i < 2; i++ {
		pn[i+4] = b*pn[i+2] - an*pn[i]
	}
	if math.Abs(pn[5]) < .0 {
		goto l35
	}
	rn = pn[4] / pn[5]
	dif = math.Abs(gin - rn)
	if dif > accurate {
		goto l34
	}
	if dif <= accurate*rn {
		goto l42
	}
l34:
	gin = rn
l35:
	for i = 0; i < 4; i++ {
		pn[i] = pn[i+2]
	}
	if math.Abs(pn[4]) < overflow {
		goto l32
	}
	for i = 0; i < 4; i++ {
		pn[i] /= overflow
	}
	goto l32
l42:
	gin = 1 - factor*gin

l50:
	return gin
}

func DiscreteGamma(alpha float64, ncat int) []float64 {
	var i int
	beta := alpha
	g := distuv.Gamma{Alpha: alpha, Beta: alpha}
	freq := make([]float64, ncat)
	r := make([]float64, ncat)
	factor := alpha / beta * float64(ncat)

	lngamma := math.Log(math.Gamma(alpha + 1))
	for i := 0; i < ncat-1; i++ {
		x := float64(i+1) / (float64(ncat))
		freq[i] = g.Quantile(x)
	}
	for i := 0; i < ncat-1; i++ {
		freq[i] = IncompleteGamma(freq[i]*beta, alpha+1, lngamma)
	}
	r[0] = freq[0] * factor
	r[ncat-1] = (1.0 - freq[ncat-2]) * factor
	for i = 1; i < ncat-1; i++ {
		r[i] = (freq[i] - freq[i-1]) * factor
	}

	return r
}

func GenerateRates(nsites int, gamma bool, alpha float64, ncat int, discrete bool) (rates []float64, categories []int) {
	rates = make([]float64, nsites)
	categories = make([]int, nsites)

	if !discrete {
		g := distuv.Gamma{Alpha: alpha, Beta: alpha}
		for i := 0; i < nsites; i++ {
			rates[i] = g.Rand()
		}
	} else if !gamma || ncat < 2 {
		for i := range rates {
			rates[i] = 1.0
			categories[i] = 0
		}
		return
	} else {
		discreteRates := DiscreteGamma(alpha, ncat)
		for i := 0; i < nsites; i++ {
			rcat := rand.Intn(ncat)
			rates[i] = discreteRates[rcat]
			categories[i] = rcat
		}
	}
	return
}
