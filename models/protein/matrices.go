package protein

import (
	"gonum.org/v1/gonum/mat"
)

// Code taken from FastME

/*********************************************************/

/* Dayhoff's model data
 * Dayhoff, M.O., Schwartz, R.M., Orcutt, B.C. (1978)
 * "A model of evolutionary change in proteins."
 * Dayhoff, M.O.(ed.) Atlas of Protein Sequence Structur., Vol5, Suppl3.
 * National Biomedical Research Foundation, Washington DC, pp.345-352. */
func DayoffMats() (dmat *mat.Dense, pi []float64) {
	var i, j, naa int
	naa = 20
	m := make([]float64, naa*naa)
	pi = make([]float64, naa)

	m[1*20+0] = 27.00
	m[2*20+0] = 98.00
	m[2*20+1] = 32.00
	m[3*20+0] = 120.00
	m[3*20+1] = 0.00
	m[3*20+2] = 905.00
	m[4*20+0] = 36.00
	m[4*20+1] = 23.00
	m[4*20+2] = 0.00
	m[4*20+3] = 0.00
	m[5*20+0] = 89.00
	m[5*20+1] = 246.00
	m[5*20+2] = 103.00
	m[5*20+3] = 134.00
	m[5*20+4] = 0.00
	m[6*20+0] = 198.00
	m[6*20+1] = 1.00
	m[6*20+2] = 148.00
	m[6*20+3] = 1153.00
	m[6*20+4] = 0.00
	m[6*20+5] = 716.00
	m[7*20+0] = 240.00
	m[7*20+1] = 9.00
	m[7*20+2] = 139.00
	m[7*20+3] = 125.00
	m[7*20+4] = 11.00
	m[7*20+5] = 28.00
	m[7*20+6] = 81.00
	m[8*20+0] = 23.00
	m[8*20+1] = 240.00
	m[8*20+2] = 535.00
	m[8*20+3] = 86.00
	m[8*20+4] = 28.00
	m[8*20+5] = 606.00
	m[8*20+6] = 43.00
	m[8*20+7] = 10.00
	m[9*20+0] = 65.00
	m[9*20+1] = 64.00
	m[9*20+2] = 77.00
	m[9*20+3] = 24.00
	m[9*20+4] = 44.00
	m[9*20+5] = 18.00
	m[9*20+6] = 61.00
	m[9*20+7] = 0.00
	m[9*20+8] = 7.00
	m[10*20+0] = 41.00
	m[10*20+1] = 15.00
	m[10*20+2] = 34.00
	m[10*20+3] = 0.00
	m[10*20+4] = 0.00
	m[10*20+5] = 73.00
	m[10*20+6] = 11.00
	m[10*20+7] = 7.00
	m[10*20+8] = 44.00
	m[10*20+9] = 257.00
	m[11*20+0] = 26.00
	m[11*20+1] = 464.00
	m[11*20+2] = 318.00
	m[11*20+3] = 71.00
	m[11*20+4] = 0.00
	m[11*20+5] = 153.00
	m[11*20+6] = 83.00
	m[11*20+7] = 27.00
	m[11*20+8] = 26.00
	m[11*20+9] = 46.00
	m[11*20+10] = 18.00
	m[12*20+0] = 72.00
	m[12*20+1] = 90.00
	m[12*20+2] = 1.00
	m[12*20+3] = 0.00
	m[12*20+4] = 0.00
	m[12*20+5] = 114.00
	m[12*20+6] = 30.00
	m[12*20+7] = 17.00
	m[12*20+8] = 0.00
	m[12*20+9] = 336.00
	m[12*20+10] = 527.00
	m[12*20+11] = 243.00
	m[13*20+0] = 18.00
	m[13*20+1] = 14.00
	m[13*20+2] = 14.00
	m[13*20+3] = 0.00
	m[13*20+4] = 0.00
	m[13*20+5] = 0.00
	m[13*20+6] = 0.00
	m[13*20+7] = 15.00
	m[13*20+8] = 48.00
	m[13*20+9] = 196.00
	m[13*20+10] = 157.00
	m[13*20+11] = 0.00
	m[13*20+12] = 92.00
	m[14*20+0] = 250.00
	m[14*20+1] = 103.00
	m[14*20+2] = 42.00
	m[14*20+3] = 13.00
	m[14*20+4] = 19.00
	m[14*20+5] = 153.00
	m[14*20+6] = 51.00
	m[14*20+7] = 34.00
	m[14*20+8] = 94.00
	m[14*20+9] = 12.00
	m[14*20+10] = 32.00
	m[14*20+11] = 33.00
	m[14*20+12] = 17.00
	m[14*20+13] = 11.00
	m[15*20+0] = 409.00
	m[15*20+1] = 154.00
	m[15*20+2] = 495.00
	m[15*20+3] = 95.00
	m[15*20+4] = 161.00
	m[15*20+5] = 56.00
	m[15*20+6] = 79.00
	m[15*20+7] = 234.00
	m[15*20+8] = 35.00
	m[15*20+9] = 24.00
	m[15*20+10] = 17.00
	m[15*20+11] = 96.00
	m[15*20+12] = 62.00
	m[15*20+13] = 46.00
	m[15*20+14] = 245.00
	m[16*20+0] = 371.00
	m[16*20+1] = 26.00
	m[16*20+2] = 229.00
	m[16*20+3] = 66.00
	m[16*20+4] = 16.00
	m[16*20+5] = 53.00
	m[16*20+6] = 34.00
	m[16*20+7] = 30.00
	m[16*20+8] = 22.00
	m[16*20+9] = 192.00
	m[16*20+10] = 33.00
	m[16*20+11] = 136.00
	m[16*20+12] = 104.00
	m[16*20+13] = 13.00
	m[16*20+14] = 78.00
	m[16*20+15] = 550.00
	m[17*20+0] = 0.00
	m[17*20+1] = 201.00
	m[17*20+2] = 23.00
	m[17*20+3] = 0.00
	m[17*20+4] = 0.00
	m[17*20+5] = 0.00
	m[17*20+6] = 0.00
	m[17*20+7] = 0.00
	m[17*20+8] = 27.00
	m[17*20+9] = 0.00
	m[17*20+10] = 46.00
	m[17*20+11] = 0.00
	m[17*20+12] = 0.00
	m[17*20+13] = 76.00
	m[17*20+14] = 0.00
	m[17*20+15] = 75.00
	m[17*20+16] = 0.00
	m[18*20+0] = 24.00
	m[18*20+1] = 8.00
	m[18*20+2] = 95.00
	m[18*20+3] = 0.00
	m[18*20+4] = 96.00
	m[18*20+5] = 0.00
	m[18*20+6] = 22.00
	m[18*20+7] = 0.00
	m[18*20+8] = 127.00
	m[18*20+9] = 37.00
	m[18*20+10] = 28.00
	m[18*20+11] = 13.00
	m[18*20+12] = 0.00
	m[18*20+13] = 698.00
	m[18*20+14] = 0.00
	m[18*20+15] = 34.00
	m[18*20+16] = 42.00
	m[18*20+17] = 61.00
	m[19*20+0] = 208.00
	m[19*20+1] = 24.00
	m[19*20+2] = 15.00
	m[19*20+3] = 18.00
	m[19*20+4] = 49.00
	m[19*20+5] = 35.00
	m[19*20+6] = 37.00
	m[19*20+7] = 54.00
	m[19*20+8] = 44.00
	m[19*20+9] = 889.00
	m[19*20+10] = 175.00
	m[19*20+11] = 10.00
	m[19*20+12] = 258.00
	m[19*20+13] = 12.00
	m[19*20+14] = 48.00
	m[19*20+15] = 30.00
	m[19*20+16] = 157.00
	m[19*20+17] = 0.00
	m[19*20+18] = 28.00

	for i = 0; i < naa; i++ {
		for j = 0; j < i; j++ {
			m[j*naa+i] = m[i*naa+j]
		}
	}

	pi[0] = 0.087127
	pi[1] = 0.040904
	pi[2] = 0.040432
	pi[3] = 0.046872
	pi[4] = 0.033474
	pi[5] = 0.038255
	pi[6] = 0.049530
	pi[7] = 0.088612
	pi[8] = 0.033618
	pi[9] = 0.036886
	pi[10] = 0.085357
	pi[11] = 0.080482
	pi[12] = 0.014753
	pi[13] = 0.039772
	pi[14] = 0.050680
	pi[15] = 0.069577
	pi[16] = 0.058542
	pi[17] = 0.010494
	pi[18] = 0.029916
	pi[19] = 0.064718

	dmat = mat.NewDense(naa, naa, m)

	return
}

/*********************************************************/

/* JTT's model data
 * D.T.Jones, W.R.Taylor and J.M.Thornton
 * "The rapid generation of mutation data matrices from protein sequences"
 * CABIOS  vol.8 no.3 1992 pp275-282 */
func JTTMats() (dmat *mat.Dense, pi []float64) {
	var i, j, naa int

	naa = 20

	m := make([]float64, naa*naa)
	pi = make([]float64, naa)

	m[1*20+0] = 58.00
	m[2*20+0] = 54.00
	m[2*20+1] = 45.00
	m[3*20+0] = 81.00
	m[3*20+1] = 16.00
	m[3*20+2] = 528.00
	m[4*20+0] = 56.00
	m[4*20+1] = 113.00
	m[4*20+2] = 34.00
	m[4*20+3] = 10.00
	m[5*20+0] = 57.00
	m[5*20+1] = 310.00
	m[5*20+2] = 86.00
	m[5*20+3] = 49.00
	m[5*20+4] = 9.00
	m[6*20+0] = 105.00
	m[6*20+1] = 29.00
	m[6*20+2] = 58.00
	m[6*20+3] = 767.00
	m[6*20+4] = 5.00
	m[6*20+5] = 323.00
	m[7*20+0] = 179.00
	m[7*20+1] = 137.00
	m[7*20+2] = 81.00
	m[7*20+3] = 130.00
	m[7*20+4] = 59.00
	m[7*20+5] = 26.00
	m[7*20+6] = 119.00
	m[8*20+0] = 27.00
	m[8*20+1] = 328.00
	m[8*20+2] = 391.00
	m[8*20+3] = 112.00
	m[8*20+4] = 69.00
	m[8*20+5] = 597.00
	m[8*20+6] = 26.00
	m[8*20+7] = 23.00
	m[9*20+0] = 36.00
	m[9*20+1] = 22.00
	m[9*20+2] = 47.00
	m[9*20+3] = 11.00
	m[9*20+4] = 17.00
	m[9*20+5] = 9.00
	m[9*20+6] = 12.00
	m[9*20+7] = 6.00
	m[9*20+8] = 16.00
	m[10*20+0] = 30.00
	m[10*20+1] = 38.00
	m[10*20+2] = 12.00
	m[10*20+3] = 7.00
	m[10*20+4] = 23.00
	m[10*20+5] = 72.00
	m[10*20+6] = 9.00
	m[10*20+7] = 6.00
	m[10*20+8] = 56.00
	m[10*20+9] = 229.00
	m[11*20+0] = 35.00
	m[11*20+1] = 646.00
	m[11*20+2] = 263.00
	m[11*20+3] = 26.00
	m[11*20+4] = 7.00
	m[11*20+5] = 292.00
	m[11*20+6] = 181.00
	m[11*20+7] = 27.00
	m[11*20+8] = 45.00
	m[11*20+9] = 21.00
	m[11*20+10] = 14.00
	m[12*20+0] = 54.00
	m[12*20+1] = 44.00
	m[12*20+2] = 30.00
	m[12*20+3] = 15.00
	m[12*20+4] = 31.00
	m[12*20+5] = 43.00
	m[12*20+6] = 18.00
	m[12*20+7] = 14.00
	m[12*20+8] = 33.00
	m[12*20+9] = 479.00
	m[12*20+10] = 388.00
	m[12*20+11] = 65.00
	m[13*20+0] = 15.00
	m[13*20+1] = 5.00
	m[13*20+2] = 10.00
	m[13*20+3] = 4.00
	m[13*20+4] = 78.00
	m[13*20+5] = 4.00
	m[13*20+6] = 5.00
	m[13*20+7] = 5.00
	m[13*20+8] = 40.00
	m[13*20+9] = 89.00
	m[13*20+10] = 248.00
	m[13*20+11] = 4.00
	m[13*20+12] = 43.00
	m[14*20+0] = 194.00
	m[14*20+1] = 74.00
	m[14*20+2] = 15.00
	m[14*20+3] = 15.00
	m[14*20+4] = 14.00
	m[14*20+5] = 164.00
	m[14*20+6] = 18.00
	m[14*20+7] = 24.00
	m[14*20+8] = 115.00
	m[14*20+9] = 10.00
	m[14*20+10] = 102.00
	m[14*20+11] = 21.00
	m[14*20+12] = 16.00
	m[14*20+13] = 17.00
	m[15*20+0] = 378.00
	m[15*20+1] = 101.00
	m[15*20+2] = 503.00
	m[15*20+3] = 59.00
	m[15*20+4] = 223.00
	m[15*20+5] = 53.00
	m[15*20+6] = 30.00
	m[15*20+7] = 201.00
	m[15*20+8] = 73.00
	m[15*20+9] = 40.00
	m[15*20+10] = 59.00
	m[15*20+11] = 47.00
	m[15*20+12] = 29.00
	m[15*20+13] = 92.00
	m[15*20+14] = 285.00
	m[16*20+0] = 475.00
	m[16*20+1] = 64.00
	m[16*20+2] = 232.00
	m[16*20+3] = 38.00
	m[16*20+4] = 42.00
	m[16*20+5] = 51.00
	m[16*20+6] = 32.00
	m[16*20+7] = 33.00
	m[16*20+8] = 46.00
	m[16*20+9] = 245.00
	m[16*20+10] = 25.00
	m[16*20+11] = 103.00
	m[16*20+12] = 226.00
	m[16*20+13] = 12.00
	m[16*20+14] = 118.00
	m[16*20+15] = 477.00
	m[17*20+0] = 9.00
	m[17*20+1] = 126.00
	m[17*20+2] = 8.00
	m[17*20+3] = 4.00
	m[17*20+4] = 115.00
	m[17*20+5] = 18.00
	m[17*20+6] = 10.00
	m[17*20+7] = 55.00
	m[17*20+8] = 8.00
	m[17*20+9] = 9.00
	m[17*20+10] = 52.00
	m[17*20+11] = 10.00
	m[17*20+12] = 24.00
	m[17*20+13] = 53.00
	m[17*20+14] = 6.00
	m[17*20+15] = 35.00
	m[17*20+16] = 12.00
	m[18*20+0] = 11.00
	m[18*20+1] = 20.00
	m[18*20+2] = 70.00
	m[18*20+3] = 46.00
	m[18*20+4] = 209.00
	m[18*20+5] = 24.00
	m[18*20+6] = 7.00
	m[18*20+7] = 8.00
	m[18*20+8] = 573.00
	m[18*20+9] = 32.00
	m[18*20+10] = 24.00
	m[18*20+11] = 8.00
	m[18*20+12] = 18.00
	m[18*20+13] = 536.00
	m[18*20+14] = 10.00
	m[18*20+15] = 63.00
	m[18*20+16] = 21.00
	m[18*20+17] = 71.00
	m[19*20+0] = 298.00
	m[19*20+1] = 17.00
	m[19*20+2] = 16.00
	m[19*20+3] = 31.00
	m[19*20+4] = 62.00
	m[19*20+5] = 20.00
	m[19*20+6] = 45.00
	m[19*20+7] = 47.00
	m[19*20+8] = 11.00
	m[19*20+9] = 961.00
	m[19*20+10] = 180.00
	m[19*20+11] = 14.00
	m[19*20+12] = 323.00
	m[19*20+13] = 62.00
	m[19*20+14] = 23.00
	m[19*20+15] = 38.00
	m[19*20+16] = 112.00
	m[19*20+17] = 25.00
	m[19*20+18] = 16.00

	for i = 0; i < naa; i++ {
		for j = 0; j < i; j++ {
			m[j*naa+i] = m[i*naa+j]
		}
	}
	pi[0] = 0.076748
	pi[1] = 0.051691
	pi[2] = 0.042645
	pi[3] = 0.051544
	pi[4] = 0.019803
	pi[5] = 0.040752
	pi[6] = 0.061830
	pi[7] = 0.073152
	pi[8] = 0.022944
	pi[9] = 0.053761
	pi[10] = 0.091904
	pi[11] = 0.058676
	pi[12] = 0.023826
	pi[13] = 0.040126
	pi[14] = 0.050901
	pi[15] = 0.068765
	pi[16] = 0.058565
	pi[17] = 0.014261
	pi[18] = 0.032102
	pi[19] = 0.066005

	dmat = mat.NewDense(naa, naa, m)
	return
}

/*********************************************************/

func MtREVMats() (dmat *mat.Dense, pi []float64) {
	var i, j, naa int

	naa = 20

	m := make([]float64, naa*naa)
	pi = make([]float64, naa)

	m[1*20+0] = 23.18
	m[2*20+0] = 26.95
	m[2*20+1] = 13.24
	m[3*20+0] = 17.67
	m[3*20+1] = 1.90
	m[3*20+2] = 794.38
	m[4*20+0] = 59.93
	m[4*20+1] = 103.33
	m[4*20+2] = 58.94
	m[4*20+3] = 1.90
	m[5*20+0] = 1.90
	m[5*20+1] = 220.99
	m[5*20+2] = 173.56
	m[5*20+3] = 55.28
	m[5*20+4] = 75.24
	m[6*20+0] = 9.77
	m[6*20+1] = 1.90
	m[6*20+2] = 63.05
	m[6*20+3] = 583.55
	m[6*20+4] = 1.90
	m[6*20+5] = 313.56
	m[7*20+0] = 120.71
	m[7*20+1] = 23.03
	m[7*20+2] = 53.30
	m[7*20+3] = 56.77
	m[7*20+4] = 30.71
	m[7*20+5] = 6.75
	m[7*20+6] = 28.28
	m[8*20+0] = 13.90
	m[8*20+1] = 165.23
	m[8*20+2] = 496.13
	m[8*20+3] = 113.99
	m[8*20+4] = 141.49
	m[8*20+5] = 582.40
	m[8*20+6] = 49.12
	m[8*20+7] = 1.90
	m[9*20+0] = 96.49
	m[9*20+1] = 1.90
	m[9*20+2] = 27.10
	m[9*20+3] = 4.34
	m[9*20+4] = 62.73
	m[9*20+5] = 8.34
	m[9*20+6] = 3.31
	m[9*20+7] = 5.98
	m[9*20+8] = 12.26
	m[10*20+0] = 25.46
	m[10*20+1] = 15.58
	m[10*20+2] = 15.16
	m[10*20+3] = 1.90
	m[10*20+4] = 25.65
	m[10*20+5] = 39.70
	m[10*20+6] = 1.90
	m[10*20+7] = 2.41
	m[10*20+8] = 11.49
	m[10*20+9] = 329.09
	m[11*20+0] = 8.36
	m[11*20+1] = 141.40
	m[11*20+2] = 608.70
	m[11*20+3] = 2.31
	m[11*20+4] = 1.90
	m[11*20+5] = 465.58
	m[11*20+6] = 313.86
	m[11*20+7] = 22.73
	m[11*20+8] = 127.67
	m[11*20+9] = 19.57
	m[11*20+10] = 14.88
	m[12*20+0] = 141.88
	m[12*20+1] = 1.90
	m[12*20+2] = 65.41
	m[12*20+3] = 1.90
	m[12*20+4] = 6.18
	m[12*20+5] = 47.37
	m[12*20+6] = 1.90
	m[12*20+7] = 1.90
	m[12*20+8] = 11.97
	m[12*20+9] = 517.98
	m[12*20+10] = 537.53
	m[12*20+11] = 91.37
	m[13*20+0] = 6.37
	m[13*20+1] = 4.69
	m[13*20+2] = 15.20
	m[13*20+3] = 4.98
	m[13*20+4] = 70.80
	m[13*20+5] = 19.11
	m[13*20+6] = 2.67
	m[13*20+7] = 1.90
	m[13*20+8] = 48.16
	m[13*20+9] = 84.67
	m[13*20+10] = 216.06
	m[13*20+11] = 6.44
	m[13*20+12] = 90.82
	m[14*20+0] = 54.31
	m[14*20+1] = 23.64
	m[14*20+2] = 73.31
	m[14*20+3] = 13.43
	m[14*20+4] = 31.26
	m[14*20+5] = 137.29
	m[14*20+6] = 12.83
	m[14*20+7] = 1.90
	m[14*20+8] = 60.97
	m[14*20+9] = 20.63
	m[14*20+10] = 40.10
	m[14*20+11] = 50.10
	m[14*20+12] = 18.84
	m[14*20+13] = 17.31
	m[15*20+0] = 387.86
	m[15*20+1] = 6.04
	m[15*20+2] = 494.39
	m[15*20+3] = 69.02
	m[15*20+4] = 277.05
	m[15*20+5] = 54.11
	m[15*20+6] = 54.71
	m[15*20+7] = 125.93
	m[15*20+8] = 77.46
	m[15*20+9] = 47.70
	m[15*20+10] = 73.61
	m[15*20+11] = 105.79
	m[15*20+12] = 111.16
	m[15*20+13] = 64.29
	m[15*20+14] = 169.90
	m[16*20+0] = 480.72
	m[16*20+1] = 2.08
	m[16*20+2] = 238.46
	m[16*20+3] = 28.01
	m[16*20+4] = 179.97
	m[16*20+5] = 94.93
	m[16*20+6] = 14.82
	m[16*20+7] = 11.17
	m[16*20+8] = 44.78
	m[16*20+9] = 368.43
	m[16*20+10] = 126.40
	m[16*20+11] = 136.33
	m[16*20+12] = 528.17
	m[16*20+13] = 33.85
	m[16*20+14] = 128.22
	m[16*20+15] = 597.21
	m[17*20+0] = 1.90
	m[17*20+1] = 21.95
	m[17*20+2] = 10.68
	m[17*20+3] = 19.86
	m[17*20+4] = 33.60
	m[17*20+5] = 1.90
	m[17*20+6] = 1.90
	m[17*20+7] = 10.92
	m[17*20+8] = 7.08
	m[17*20+9] = 1.90
	m[17*20+10] = 32.44
	m[17*20+11] = 24.00
	m[17*20+12] = 21.71
	m[17*20+13] = 7.84
	m[17*20+14] = 4.21
	m[17*20+15] = 38.58
	m[17*20+16] = 9.99
	m[18*20+0] = 6.48
	m[18*20+1] = 1.90
	m[18*20+2] = 191.36
	m[18*20+3] = 21.21
	m[18*20+4] = 254.77
	m[18*20+5] = 38.82
	m[18*20+6] = 13.12
	m[18*20+7] = 3.21
	m[18*20+8] = 670.14
	m[18*20+9] = 25.01
	m[18*20+10] = 44.15
	m[18*20+11] = 51.17
	m[18*20+12] = 39.96
	m[18*20+13] = 465.58
	m[18*20+14] = 16.21
	m[18*20+15] = 64.92
	m[18*20+16] = 38.73
	m[18*20+17] = 26.25
	m[19*20+0] = 195.06
	m[19*20+1] = 7.64
	m[19*20+2] = 1.90
	m[19*20+3] = 1.90
	m[19*20+4] = 1.90
	m[19*20+5] = 19.00
	m[19*20+6] = 21.14
	m[19*20+7] = 2.53
	m[19*20+8] = 1.90
	m[19*20+9] = 1222.94
	m[19*20+10] = 91.67
	m[19*20+11] = 1.90
	m[19*20+12] = 387.54
	m[19*20+13] = 6.35
	m[19*20+14] = 8.23
	m[19*20+15] = 1.90
	m[19*20+16] = 204.54
	m[19*20+17] = 5.37
	m[19*20+18] = 1.90

	for i = 0; i < naa; i++ {
		for j = 0; j < i; j++ {
			m[j*naa+i] = m[i*naa+j]
		}
	}
	pi[0] = 0.072000
	pi[1] = 0.019000
	pi[2] = 0.039000
	pi[3] = 0.019000
	pi[4] = 0.006000
	pi[5] = 0.025000
	pi[6] = 0.024000
	pi[7] = 0.056000
	pi[8] = 0.028000
	pi[9] = 0.088000
	pi[10] = 0.169000
	pi[11] = 0.023000
	pi[12] = 0.054000
	pi[13] = 0.061000
	pi[14] = 0.054000
	pi[15] = 0.072000
	pi[16] = 0.086000
	pi[17] = 0.029000
	pi[18] = 0.033000
	pi[19] = 0.043000

	dmat = mat.NewDense(naa, naa, m)
	return
}

/*********************************************************/

/* LG model
 * Si Quang LE & Olivier Gascuel
 * "An improved general amino-acid replacement matrix"
 * Mol Biol Evol. 2008 Jul;25(7):1307-20. */
func LGMats() (dmat *mat.Dense, pi []float64) {
	var i, j, naa int

	naa = 20

	m := make([]float64, naa*naa)
	pi = make([]float64, naa)

	m[1*20+0] = 0.425093
	m[2*20+0] = 0.276818
	m[2*20+1] = 0.751878
	m[3*20+0] = 0.395144
	m[3*20+1] = 0.123954
	m[3*20+2] = 5.076149
	m[4*20+0] = 2.489084
	m[4*20+1] = 0.534551
	m[4*20+2] = 0.528768
	m[4*20+3] = 0.062556
	m[5*20+0] = 0.969894
	m[5*20+1] = 2.807908
	m[5*20+2] = 1.695752
	m[5*20+3] = 0.523386
	m[5*20+4] = 0.084808
	m[6*20+0] = 1.038545
	m[6*20+1] = 0.363970
	m[6*20+2] = 0.541712
	m[6*20+3] = 5.243870
	m[6*20+4] = 0.003499
	m[6*20+5] = 4.128591
	m[7*20+0] = 2.066040
	m[7*20+1] = 0.390192
	m[7*20+2] = 1.437645
	m[7*20+3] = 0.844926
	m[7*20+4] = 0.569265
	m[7*20+5] = 0.267959
	m[7*20+6] = 0.348847
	m[8*20+0] = 0.358858
	m[8*20+1] = 2.426601
	m[8*20+2] = 4.509238
	m[8*20+3] = 0.927114
	m[8*20+4] = 0.640543
	m[8*20+5] = 4.813505
	m[8*20+6] = 0.423881
	m[8*20+7] = 0.311484
	m[9*20+0] = 0.149830
	m[9*20+1] = 0.126991
	m[9*20+2] = 0.191503
	m[9*20+3] = 0.010690
	m[9*20+4] = 0.320627
	m[9*20+5] = 0.072854
	m[9*20+6] = 0.044265
	m[9*20+7] = 0.008705
	m[9*20+8] = 0.108882
	m[10*20+0] = 0.395337
	m[10*20+1] = 0.301848
	m[10*20+2] = 0.068427
	m[10*20+3] = 0.015076
	m[10*20+4] = 0.594007
	m[10*20+5] = 0.582457
	m[10*20+6] = 0.069673
	m[10*20+7] = 0.044261
	m[10*20+8] = 0.366317
	m[10*20+9] = 4.145067
	m[11*20+0] = 0.536518
	m[11*20+1] = 6.326067
	m[11*20+2] = 2.145078
	m[11*20+3] = 0.282959
	m[11*20+4] = 0.013266
	m[11*20+5] = 3.234294
	m[11*20+6] = 1.807177
	m[11*20+7] = 0.296636
	m[11*20+8] = 0.697264
	m[11*20+9] = 0.159069
	m[11*20+10] = 0.137500
	m[12*20+0] = 1.124035
	m[12*20+1] = 0.484133
	m[12*20+2] = 0.371004
	m[12*20+3] = 0.025548
	m[12*20+4] = 0.893680
	m[12*20+5] = 1.672569
	m[12*20+6] = 0.173735
	m[12*20+7] = 0.139538
	m[12*20+8] = 0.442472
	m[12*20+9] = 4.273607
	m[12*20+10] = 6.312358
	m[12*20+11] = 0.656604
	m[13*20+0] = 0.253701
	m[13*20+1] = 0.052722
	m[13*20+2] = 0.089525
	m[13*20+3] = 0.017416
	m[13*20+4] = 1.105251
	m[13*20+5] = 0.035855
	m[13*20+6] = 0.018811
	m[13*20+7] = 0.089586
	m[13*20+8] = 0.682139
	m[13*20+9] = 1.112727
	m[13*20+10] = 2.592692
	m[13*20+11] = 0.023918
	m[13*20+12] = 1.798853
	m[14*20+0] = 1.177651
	m[14*20+1] = 0.332533
	m[14*20+2] = 0.161787
	m[14*20+3] = 0.394456
	m[14*20+4] = 0.075382
	m[14*20+5] = 0.624294
	m[14*20+6] = 0.419409
	m[14*20+7] = 0.196961
	m[14*20+8] = 0.508851
	m[14*20+9] = 0.078281
	m[14*20+10] = 0.249060
	m[14*20+11] = 0.390322
	m[14*20+12] = 0.099849
	m[14*20+13] = 0.094464
	m[15*20+0] = 4.727182
	m[15*20+1] = 0.858151
	m[15*20+2] = 4.008358
	m[15*20+3] = 1.240275
	m[15*20+4] = 2.784478
	m[15*20+5] = 1.223828
	m[15*20+6] = 0.611973
	m[15*20+7] = 1.739990
	m[15*20+8] = 0.990012
	m[15*20+9] = 0.064105
	m[15*20+10] = 0.182287
	m[15*20+11] = 0.748683
	m[15*20+12] = 0.346960
	m[15*20+13] = 0.361819
	m[15*20+14] = 1.338132
	m[16*20+0] = 2.139501
	m[16*20+1] = 0.578987
	m[16*20+2] = 2.000679
	m[16*20+3] = 0.425860
	m[16*20+4] = 1.143480
	m[16*20+5] = 1.080136
	m[16*20+6] = 0.604545
	m[16*20+7] = 0.129836
	m[16*20+8] = 0.584262
	m[16*20+9] = 1.033739
	m[16*20+10] = 0.302936
	m[16*20+11] = 1.136863
	m[16*20+12] = 2.020366
	m[16*20+13] = 0.165001
	m[16*20+14] = 0.571468
	m[16*20+15] = 6.472279
	m[17*20+0] = 0.180717
	m[17*20+1] = 0.593607
	m[17*20+2] = 0.045376
	m[17*20+3] = 0.029890
	m[17*20+4] = 0.670128
	m[17*20+5] = 0.236199
	m[17*20+6] = 0.077852
	m[17*20+7] = 0.268491
	m[17*20+8] = 0.597054
	m[17*20+9] = 0.111660
	m[17*20+10] = 0.619632
	m[17*20+11] = 0.049906
	m[17*20+12] = 0.696175
	m[17*20+13] = 2.457121
	m[17*20+14] = 0.095131
	m[17*20+15] = 0.248862
	m[17*20+16] = 0.140825
	m[18*20+0] = 0.218959
	m[18*20+1] = 0.314440
	m[18*20+2] = 0.612025
	m[18*20+3] = 0.135107
	m[18*20+4] = 1.165532
	m[18*20+5] = 0.257336
	m[18*20+6] = 0.120037
	m[18*20+7] = 0.054679
	m[18*20+8] = 5.306834
	m[18*20+9] = 0.232523
	m[18*20+10] = 0.299648
	m[18*20+11] = 0.131932
	m[18*20+12] = 0.481306
	m[18*20+13] = 7.803902
	m[18*20+14] = 0.089613
	m[18*20+15] = 0.400547
	m[18*20+16] = 0.245841
	m[18*20+17] = 3.151815
	m[19*20+0] = 2.547870
	m[19*20+1] = 0.170887
	m[19*20+2] = 0.083688
	m[19*20+3] = 0.037967
	m[19*20+4] = 1.959291
	m[19*20+5] = 0.210332
	m[19*20+6] = 0.245034
	m[19*20+7] = 0.076701
	m[19*20+8] = 0.119013
	m[19*20+9] = 10.649107
	m[19*20+10] = 1.702745
	m[19*20+11] = 0.185202
	m[19*20+12] = 1.898718
	m[19*20+13] = 0.654683
	m[19*20+14] = 0.296501
	m[19*20+15] = 0.098369
	m[19*20+16] = 2.188158
	m[19*20+17] = 0.189510
	m[19*20+18] = 0.249313

	for i = 0; i < naa; i++ {
		for j = 0; j < i; j++ {
			m[j*naa+i] = m[i*naa+j]
		}
	}
	pi[0] = 0.079066
	pi[1] = 0.055941
	pi[2] = 0.041977
	pi[3] = 0.053052
	pi[4] = 0.012937
	pi[5] = 0.040767
	pi[6] = 0.071586
	pi[7] = 0.057337
	pi[8] = 0.022355
	pi[9] = 0.062157
	pi[10] = 0.099081
	pi[11] = 0.064600
	pi[12] = 0.022951
	pi[13] = 0.042302
	pi[14] = 0.044040
	pi[15] = 0.061197
	pi[16] = 0.053287
	pi[17] = 0.012066
	pi[18] = 0.034155
	pi[19] = 0.069147

	dmat = mat.NewDense(naa, naa, m)
	return
}

/*********************************************************/

/* WAG's model data
 * Simon Whelan and Nick Goldman
 * "A general empirical model of protein evolution derived from multiple
 * protein families using a maximum-likelihood approach"
 * MBE (2001) 18:691-699 */
func WAGMats() (dmat *mat.Dense, pi []float64) {
	var i, j, naa int

	naa = 20

	m := make([]float64, naa*naa)
	pi = make([]float64, naa)

	m[1*20+0] = 55.15710
	m[2*20+0] = 50.98480
	m[2*20+1] = 63.53460
	m[3*20+0] = 73.89980
	m[3*20+1] = 14.73040
	m[3*20+2] = 542.94200
	m[4*20+0] = 102.70400
	m[4*20+1] = 52.81910
	m[4*20+2] = 26.52560
	m[4*20+3] = 3.02949
	m[5*20+0] = 90.85980
	m[5*20+1] = 303.55000
	m[5*20+2] = 154.36400
	m[5*20+3] = 61.67830
	m[5*20+4] = 9.88179
	m[6*20+0] = 158.28500
	m[6*20+1] = 43.91570
	m[6*20+2] = 94.71980
	m[6*20+3] = 617.41600
	m[6*20+4] = 2.13520
	m[6*20+5] = 546.94700
	m[7*20+0] = 141.67200
	m[7*20+1] = 58.46650
	m[7*20+2] = 112.55600
	m[7*20+3] = 86.55840
	m[7*20+4] = 30.66740
	m[7*20+5] = 33.00520
	m[7*20+6] = 56.77170
	m[8*20+0] = 31.69540
	m[8*20+1] = 213.71500
	m[8*20+2] = 395.62900
	m[8*20+3] = 93.06760
	m[8*20+4] = 24.89720
	m[8*20+5] = 429.41100
	m[8*20+6] = 57.00250
	m[8*20+7] = 24.94100
	m[9*20+0] = 19.33350
	m[9*20+1] = 18.69790
	m[9*20+2] = 55.42360
	m[9*20+3] = 3.94370
	m[9*20+4] = 17.01350
	m[9*20+5] = 11.39170
	m[9*20+6] = 12.73950
	m[9*20+7] = 3.04501
	m[9*20+8] = 13.81900
	m[10*20+0] = 39.79150
	m[10*20+1] = 49.76710
	m[10*20+2] = 13.15280
	m[10*20+3] = 8.48047
	m[10*20+4] = 38.42870
	m[10*20+5] = 86.94890
	m[10*20+6] = 15.42630
	m[10*20+7] = 6.13037
	m[10*20+8] = 49.94620
	m[10*20+9] = 317.09700
	m[11*20+0] = 90.62650
	m[11*20+1] = 535.14200
	m[11*20+2] = 301.20100
	m[11*20+3] = 47.98550
	m[11*20+4] = 7.40339
	m[11*20+5] = 389.49000
	m[11*20+6] = 258.44300
	m[11*20+7] = 37.35580
	m[11*20+8] = 89.04320
	m[11*20+9] = 32.38320
	m[11*20+10] = 25.75550
	m[12*20+0] = 89.34960
	m[12*20+1] = 68.31620
	m[12*20+2] = 19.82210
	m[12*20+3] = 10.37540
	m[12*20+4] = 39.04820
	m[12*20+5] = 154.52600
	m[12*20+6] = 31.51240
	m[12*20+7] = 17.41000
	m[12*20+8] = 40.41410
	m[12*20+9] = 425.74600
	m[12*20+10] = 485.40200
	m[12*20+11] = 93.42760
	m[13*20+0] = 21.04940
	m[13*20+1] = 10.27110
	m[13*20+2] = 9.61621
	m[13*20+3] = 4.67304
	m[13*20+4] = 39.80200
	m[13*20+5] = 9.99208
	m[13*20+6] = 8.11339
	m[13*20+7] = 4.99310
	m[13*20+8] = 67.93710
	m[13*20+9] = 105.94700
	m[13*20+10] = 211.51700
	m[13*20+11] = 8.88360
	m[13*20+12] = 119.06300
	m[14*20+0] = 143.85500
	m[14*20+1] = 67.94890
	m[14*20+2] = 19.50810
	m[14*20+3] = 42.39840
	m[14*20+4] = 10.94040
	m[14*20+5] = 93.33720
	m[14*20+6] = 68.23550
	m[14*20+7] = 24.35700
	m[14*20+8] = 69.61980
	m[14*20+9] = 9.99288
	m[14*20+10] = 41.58440
	m[14*20+11] = 55.68960
	m[14*20+12] = 17.13290
	m[14*20+13] = 16.14440
	m[15*20+0] = 337.07900
	m[15*20+1] = 122.41900
	m[15*20+2] = 397.42300
	m[15*20+3] = 107.17600
	m[15*20+4] = 140.76600
	m[15*20+5] = 102.88700
	m[15*20+6] = 70.49390
	m[15*20+7] = 134.18200
	m[15*20+8] = 74.01690
	m[15*20+9] = 31.94400
	m[15*20+10] = 34.47390
	m[15*20+11] = 96.71300
	m[15*20+12] = 49.39050
	m[15*20+13] = 54.59310
	m[15*20+14] = 161.32800
	m[16*20+0] = 212.11100
	m[16*20+1] = 55.44130
	m[16*20+2] = 203.00600
	m[16*20+3] = 37.48660
	m[16*20+4] = 51.29840
	m[16*20+5] = 85.79280
	m[16*20+6] = 82.27650
	m[16*20+7] = 22.58330
	m[16*20+8] = 47.33070
	m[16*20+9] = 145.81600
	m[16*20+10] = 32.66220
	m[16*20+11] = 138.69800
	m[16*20+12] = 151.61200
	m[16*20+13] = 17.19030
	m[16*20+14] = 79.53840
	m[16*20+15] = 437.80200
	m[17*20+0] = 11.31330
	m[17*20+1] = 116.39200
	m[17*20+2] = 7.19167
	m[17*20+3] = 12.97670
	m[17*20+4] = 71.70700
	m[17*20+5] = 21.57370
	m[17*20+6] = 15.65570
	m[17*20+7] = 33.69830
	m[17*20+8] = 26.25690
	m[17*20+9] = 21.24830
	m[17*20+10] = 66.53090
	m[17*20+11] = 13.75050
	m[17*20+12] = 51.57060
	m[17*20+13] = 152.96400
	m[17*20+14] = 13.94050
	m[17*20+15] = 52.37420
	m[17*20+16] = 11.08640
	m[18*20+0] = 24.07350
	m[18*20+1] = 38.15330
	m[18*20+2] = 108.60000
	m[18*20+3] = 32.57110
	m[18*20+4] = 54.38330
	m[18*20+5] = 22.77100
	m[18*20+6] = 19.63030
	m[18*20+7] = 10.36040
	m[18*20+8] = 387.34400
	m[18*20+9] = 42.01700
	m[18*20+10] = 39.86180
	m[18*20+11] = 13.32640
	m[18*20+12] = 42.84370
	m[18*20+13] = 645.42800
	m[18*20+14] = 21.60460
	m[18*20+15] = 78.69930
	m[18*20+16] = 29.11480
	m[18*20+17] = 248.53900
	m[19*20+0] = 200.60100
	m[19*20+1] = 25.18490
	m[19*20+2] = 19.62460
	m[19*20+3] = 15.23350
	m[19*20+4] = 100.21400
	m[19*20+5] = 30.12810
	m[19*20+6] = 58.87310
	m[19*20+7] = 18.72470
	m[19*20+8] = 11.83580
	m[19*20+9] = 782.13000
	m[19*20+10] = 180.03400
	m[19*20+11] = 30.54340
	m[19*20+12] = 205.84500
	m[19*20+13] = 64.98920
	m[19*20+14] = 31.48870
	m[19*20+15] = 23.27390
	m[19*20+16] = 138.82300
	m[19*20+17] = 36.53690
	m[19*20+18] = 31.47300

	for i = 0; i < naa; i++ {
		for j = 0; j < i; j++ {
			m[j*naa+i] = m[i*naa+j]
		}
	}
	pi[0] = 0.0866279
	pi[1] = 0.043972
	pi[2] = 0.0390894
	pi[3] = 0.0570451
	pi[4] = 0.0193078
	pi[5] = 0.0367281
	pi[6] = 0.0580589
	pi[7] = 0.0832518
	pi[8] = 0.0244313
	pi[9] = 0.048466
	pi[10] = 0.086209
	pi[11] = 0.0620286
	pi[12] = 0.0195027
	pi[13] = 0.0384319
	pi[14] = 0.0457631
	pi[15] = 0.0695179
	pi[16] = 0.0610127
	pi[17] = 0.0143859
	pi[18] = 0.0352742
	pi[19] = 0.0708956

	dmat = mat.NewDense(naa, naa, m)
	return
}

// From PhyML code
// Nickle DC, Heath L, Jensen MA, Gilbert PB, Mullins JI, Kosakovsky Pond SL.
// HIV-Specific Probabilistic Models of Protein Evolution.
// PLoS ONE. 2007 Jun 6;2:e503.
// [thanks to Sergei L. Kosakovsky]
// Translated from HYPHY to Phyml format by Federico Abascal.
func HIVBMats() (dmat *mat.Dense, pi []float64) {
	var i, j, naa int
	naa = 20
	m := make([]float64, naa*naa)
	pi = make([]float64, naa)

	m[1*20+0] = 0.307507
	m[2*20+0] = 0.005
	m[2*20+1] = 0.295543
	m[3*20+0] = 1.45504
	m[3*20+1] = 0.005
	m[3*20+2] = 17.6612
	m[4*20+0] = 0.123758
	m[4*20+1] = 0.351721
	m[4*20+2] = 0.0860642
	m[4*20+3] = 0.005
	m[5*20+0] = 0.0551128
	m[5*20+1] = 3.4215
	m[5*20+2] = 0.672052
	m[5*20+3] = 0.005
	m[5*20+4] = 0.005
	m[6*20+0] = 1.48135
	m[6*20+1] = 0.0749218
	m[6*20+2] = 0.0792633
	m[6*20+3] = 10.5872
	m[6*20+4] = 0.005
	m[6*20+5] = 2.5602
	m[7*20+0] = 2.13536
	m[7*20+1] = 3.65345
	m[7*20+2] = 0.323401
	m[7*20+3] = 2.83806
	m[7*20+4] = 0.897871
	m[7*20+5] = 0.0619137
	m[7*20+6] = 3.92775
	m[8*20+0] = 0.0847613
	m[8*20+1] = 9.04044
	m[8*20+2] = 7.64585
	m[8*20+3] = 1.9169
	m[8*20+4] = 0.240073
	m[8*20+5] = 7.05545
	m[8*20+6] = 0.11974
	m[8*20+7] = 0.005
	m[9*20+0] = 0.005
	m[9*20+1] = 0.677289
	m[9*20+2] = 0.680565
	m[9*20+3] = 0.0176792
	m[9*20+4] = 0.005
	m[9*20+5] = 0.005
	m[9*20+6] = 0.00609079
	m[9*20+7] = 0.005
	m[9*20+8] = 0.103111
	m[10*20+0] = 0.215256
	m[10*20+1] = 0.701427
	m[10*20+2] = 0.005
	m[10*20+3] = 0.00876048
	m[10*20+4] = 0.129777
	m[10*20+5] = 1.49456
	m[10*20+6] = 0.005
	m[10*20+7] = 0.005
	m[10*20+8] = 1.74171
	m[10*20+9] = 5.95879
	m[11*20+0] = 0.005
	m[11*20+1] = 20.45
	m[11*20+2] = 7.90443
	m[11*20+3] = 0.005
	m[11*20+4] = 0.005
	m[11*20+5] = 6.54737
	m[11*20+6] = 4.61482
	m[11*20+7] = 0.521705
	m[11*20+8] = 0.005
	m[11*20+9] = 0.322319
	m[11*20+10] = 0.0814995
	m[12*20+0] = 0.0186643
	m[12*20+1] = 2.51394
	m[12*20+2] = 0.005
	m[12*20+3] = 0.005
	m[12*20+4] = 0.005
	m[12*20+5] = 0.303676
	m[12*20+6] = 0.175789
	m[12*20+7] = 0.005
	m[12*20+8] = 0.005
	m[12*20+9] = 11.2065
	m[12*20+10] = 5.31961
	m[12*20+11] = 1.28246
	m[13*20+0] = 0.0141269
	m[13*20+1] = 0.005
	m[13*20+2] = 0.005
	m[13*20+3] = 0.005
	m[13*20+4] = 9.29815
	m[13*20+5] = 0.005
	m[13*20+6] = 0.005
	m[13*20+7] = 0.291561
	m[13*20+8] = 0.145558
	m[13*20+9] = 3.39836
	m[13*20+10] = 8.52484
	m[13*20+11] = 0.0342658
	m[13*20+12] = 0.188025
	m[14*20+0] = 2.12217
	m[14*20+1] = 1.28355
	m[14*20+2] = 0.00739578
	m[14*20+3] = 0.0342658
	m[14*20+4] = 0.005
	m[14*20+5] = 4.47211
	m[14*20+6] = 0.0120226
	m[14*20+7] = 0.005
	m[14*20+8] = 2.45318
	m[14*20+9] = 0.0410593
	m[14*20+10] = 2.07757
	m[14*20+11] = 0.0313862
	m[14*20+12] = 0.005
	m[14*20+13] = 0.005
	m[15*20+0] = 2.46633
	m[15*20+1] = 3.4791
	m[15*20+2] = 13.1447
	m[15*20+3] = 0.52823
	m[15*20+4] = 4.69314
	m[15*20+5] = 0.116311
	m[15*20+6] = 0.005
	m[15*20+7] = 4.38041
	m[15*20+8] = 0.382747
	m[15*20+9] = 1.21803
	m[15*20+10] = 0.927656
	m[15*20+11] = 0.504111
	m[15*20+12] = 0.005
	m[15*20+13] = 0.956472
	m[15*20+14] = 5.37762
	m[16*20+0] = 15.9183
	m[16*20+1] = 2.86868
	m[16*20+2] = 6.88667
	m[16*20+3] = 0.274724
	m[16*20+4] = 0.739969
	m[16*20+5] = 0.243589
	m[16*20+6] = 0.289774
	m[16*20+7] = 0.369615
	m[16*20+8] = 0.711594
	m[16*20+9] = 8.61217
	m[16*20+10] = 0.0437673
	m[16*20+11] = 4.67142
	m[16*20+12] = 4.94026
	m[16*20+13] = 0.0141269
	m[16*20+14] = 2.01417
	m[16*20+15] = 8.93107
	m[17*20+0] = 0.005
	m[17*20+1] = 0.991338
	m[17*20+2] = 0.005
	m[17*20+3] = 0.005
	m[17*20+4] = 2.63277
	m[17*20+5] = 0.026656
	m[17*20+6] = 0.005
	m[17*20+7] = 1.21674
	m[17*20+8] = 0.0695179
	m[17*20+9] = 0.005
	m[17*20+10] = 0.748843
	m[17*20+11] = 0.005
	m[17*20+12] = 0.089078
	m[17*20+13] = 0.829343
	m[17*20+14] = 0.0444506
	m[17*20+15] = 0.0248728
	m[17*20+16] = 0.005
	m[18*20+0] = 0.005
	m[18*20+1] = 0.00991826
	m[18*20+2] = 1.76417
	m[18*20+3] = 0.674653
	m[18*20+4] = 7.57932
	m[18*20+5] = 0.113033
	m[18*20+6] = 0.0792633
	m[18*20+7] = 0.005
	m[18*20+8] = 18.6943
	m[18*20+9] = 0.148168
	m[18*20+10] = 0.111986
	m[18*20+11] = 0.005
	m[18*20+12] = 0.005
	m[18*20+13] = 15.34
	m[18*20+14] = 0.0304381
	m[18*20+15] = 0.648024
	m[18*20+16] = 0.105652
	m[18*20+17] = 1.28022
	m[19*20+0] = 7.61428
	m[19*20+1] = 0.0812454
	m[19*20+2] = 0.026656
	m[19*20+3] = 1.04793
	m[19*20+4] = 0.420027
	m[19*20+5] = 0.0209153
	m[19*20+6] = 1.02847
	m[19*20+7] = 0.953155
	m[19*20+8] = 0.005
	m[19*20+9] = 17.7389
	m[19*20+10] = 1.41036
	m[19*20+11] = 0.265829
	m[19*20+12] = 6.8532
	m[19*20+13] = 0.723274
	m[19*20+14] = 0.005
	m[19*20+15] = 0.0749218
	m[19*20+16] = 0.709226
	m[19*20+17] = 0.005
	m[19*20+18] = 0.0410593

	for i = 0; i < naa; i++ {
		for j = 0; j < i; j++ {
			m[j*naa+i] = m[i*naa+j]
		}
	}

	pi[0] = 0.060490222
	pi[1] = 0.066039665
	pi[2] = 0.044127815
	pi[3] = 0.042109048
	pi[4] = 0.020075899
	pi[5] = 0.053606488
	pi[6] = 0.071567447
	pi[7] = 0.072308239
	pi[8] = 0.022293943
	pi[9] = 0.069730629
	pi[10] = 0.098851122
	pi[11] = 0.056968211
	pi[12] = 0.019768318
	pi[13] = 0.028809447
	pi[14] = 0.046025282
	pi[15] = 0.05060433
	pi[16] = 0.053636813
	pi[17] = 0.033011601
	pi[18] = 0.028350243
	pi[19] = 0.061625237

	dmat = mat.NewDense(naa, naa, m)
	return
}

// Antibody-Specific Model of Amino Acid Substitution
// From AB Model publication (supplementary material online)
// Alexander Mirsky, Linda Kazandjian, Maria Anisimova,
// Antibody-Specific Model of Amino Acid Substitution for Immunological Inferences from Alignments of Antibody Sequences,
// Molecular Biology and Evolution, Volume 32, Issue 3, March 2015, Pages 806–819,
// https://doi.org/10.1093/molbev/msu340
func ABMats() (dmat *mat.Dense, pi []float64) {
	var i, j, naa int
	naa = 20
	m := make([]float64, naa*naa)
	pi = make([]float64, naa)

	m[1*20+0] = 1.784266e-01
	m[2*20+0] = 9.291290e-02
	m[2*20+1] = 7.829130e-01
	m[3*20+0] = 1.241095e+00
	m[3*20+1] = 5.795374e-02
	m[3*20+2] = 7.185182e+00
	m[4*20+0] = 8.929181e-03
	m[4*20+1] = 1.821885e-01
	m[4*20+2] = 1.374268e-06
	m[4*20+3] = 2.340019e-02
	m[5*20+0] = 1.992269e-01
	m[5*20+1] = 1.923901e+00
	m[5*20+2] = 8.705989e-02
	m[5*20+3] = 1.843856e-01
	m[5*20+4] = 1.046446e-08
	m[6*20+0] = 9.521821e-01
	m[6*20+1] = 6.273863e-02
	m[6*20+2] = 5.038373e-01
	m[6*20+3] = 7.426619e+00
	m[6*20+4] = 7.519215e-11
	m[6*20+5] = 3.691671e+00
	m[7*20+0] = 1.851951e+00
	m[7*20+1] = 1.089400e+00
	m[7*20+2] = 4.868901e-01
	m[7*20+3] = 2.112400e+00
	m[7*20+4] = 5.891123e-02
	m[7*20+5] = 5.516340e-02
	m[7*20+6] = 1.389370e+00
	m[8*20+0] = 5.241316e+00
	m[8*20+1] = 1.049550e+01
	m[8*20+2] = 1.405444e+01
	m[8*20+3] = 1.126995e+01
	m[8*20+4] = 3.963388e+00
	m[8*20+5] = 8.908434e+00
	m[8*20+6] = 7.298080e+00
	m[8*20+7] = 9.139518e+00
	m[9*20+0] = 1.140412e-01
	m[9*20+1] = 3.245175e-01
	m[9*20+2] = 1.762721e+00
	m[9*20+3] = 3.916999e-02
	m[9*20+4] = 6.594967e-04
	m[9*20+5] = 6.712736e-06
	m[9*20+6] = 1.029959e-04
	m[9*20+7] = 3.560482e-02
	m[9*20+8] = 4.706586e+00
	m[10*20+0] = 6.969101e-02
	m[10*20+1] = 3.932002e-01
	m[10*20+2] = 2.769442e-02
	m[10*20+3] = 3.020502e-02
	m[10*20+4] = 6.079219e-03
	m[10*20+5] = 6.802781e-01
	m[10*20+6] = 1.283121e-03
	m[10*20+7] = 2.157936e-02
	m[10*20+8] = 5.879103e+00
	m[10*20+9] = 1.601123e+00
	m[11*20+0] = 7.388355e-02
	m[11*20+1] = 7.549240e+00
	m[11*20+2] = 6.190318e+00
	m[11*20+3] = 6.622772e-02
	m[11*20+4] = 3.722878e-16
	m[11*20+5] = 3.030805e+00
	m[11*20+6] = 3.608816e+00
	m[11*20+7] = 5.504400e-02
	m[11*20+8] = 1.455741e+00
	m[11*20+9] = 5.059793e-01
	m[11*20+10] = 2.158451e-02
	m[12*20+0] = 6.299271e-02
	m[12*20+1] = 3.362326e-01
	m[12*20+2] = 3.972173e-02
	m[12*20+3] = 3.357577e-02
	m[12*20+4] = 7.213178e-03
	m[12*20+5] = 1.233336e-03
	m[12*20+6] = 7.659566e-02
	m[12*20+7] = 2.187264e-02
	m[12*20+8] = 2.298295e+00
	m[12*20+9] = 1.096748e+01
	m[12*20+10] = 5.647985e+00
	m[12*20+11] = 1.238634e+00
	m[13*20+0] = 1.130146e-01
	m[13*20+1] = 8.208677e-02
	m[13*20+2] = 1.955446e-01
	m[13*20+3] = 1.031734e-01
	m[13*20+4] = 1.993818e-01
	m[13*20+5] = 1.496610e-03
	m[13*20+6] = 5.288625e-02
	m[13*20+7] = 1.984772e-01
	m[13*20+8] = 5.642309e+00
	m[13*20+9] = 2.714705e+00
	m[13*20+10] = 3.390618e+00
	m[13*20+11] = 4.649035e-03
	m[13*20+12] = 3.947940e+00
	m[14*20+0] = 1.800713e+00
	m[14*20+1] = 3.498713e-01
	m[14*20+2] = 7.342554e-03
	m[14*20+3] = 1.509482e-01
	m[14*20+4] = 4.878395e-03
	m[14*20+5] = 7.426909e-01
	m[14*20+6] = 2.889815e-02
	m[14*20+7] = 7.915056e-02
	m[14*20+8] = 1.049496e+01
	m[14*20+9] = 5.016568e-02
	m[14*20+10] = 1.149931e+00
	m[14*20+11] = 9.948994e-03
	m[14*20+12] = 7.417279e-02
	m[14*20+13] = 3.556198e-01
	m[15*20+0] = 9.988358e-01
	m[15*20+1] = 1.926435e+00
	m[15*20+2] = 7.348346e+00
	m[15*20+3] = 5.822988e-01
	m[15*20+4] = 2.639482e-01
	m[15*20+5] = 5.906405e-04
	m[15*20+6] = 6.776709e-02
	m[15*20+7] = 9.984215e-01
	m[15*20+8] = 5.439116e+00
	m[15*20+9] = 6.007607e-01
	m[15*20+10] = 1.580539e-01
	m[15*20+11] = 8.688405e-02
	m[15*20+12] = 1.861354e-02
	m[15*20+13] = 9.813064e-01
	m[15*20+14] = 1.284651e+00
	m[16*20+0] = 2.912317e+00
	m[16*20+1] = 1.135258e+00
	m[16*20+2] = 2.147175e+00
	m[16*20+3] = 1.516881e-01
	m[16*20+4] = 3.225214e-06
	m[16*20+5] = 1.202094e-01
	m[16*20+6] = 6.016624e-02
	m[16*20+7] = 7.862767e-02
	m[16*20+8] = 3.443285e+00
	m[16*20+9] = 3.087152e+00
	m[16*20+10] = 5.702792e-01
	m[16*20+11] = 1.039298e+00
	m[16*20+12] = 1.415612e+00
	m[16*20+13] = 3.674486e-02
	m[16*20+14] = 9.057112e-01
	m[16*20+15] = 3.058575e+00
	m[17*20+0] = 7.939549e-02
	m[17*20+1] = 5.724286e-01
	m[17*20+2] = 7.310937e-04
	m[17*20+3] = 1.423897e-02
	m[17*20+4] = 4.440833e-01
	m[17*20+5] = 4.332983e-05
	m[17*20+6] = 2.252612e-02
	m[17*20+7] = 1.386853e-01
	m[17*20+8] = 7.013890e+00
	m[17*20+9] = 6.318748e-02
	m[17*20+10] = 3.378544e-01
	m[17*20+11] = 8.024263e-03
	m[17*20+12] = 1.011149e-01
	m[17*20+13] = 2.199856e-01
	m[17*20+14] = 5.516074e-03
	m[17*20+15] = 1.385142e-01
	m[17*20+16] = 1.412361e-02
	m[18*20+0] = 1.433528e-01
	m[18*20+1] = 1.711315e-01
	m[18*20+2] = 2.622763e+00
	m[18*20+3] = 9.078338e-01
	m[18*20+4] = 7.741612e-01
	m[18*20+5] = 2.737091e-02
	m[18*20+6] = 1.240642e-01
	m[18*20+7] = 2.295842e-01
	m[18*20+8] = 2.055414e+01
	m[18*20+9] = 2.903165e-01
	m[18*20+10] = 1.521320e-01
	m[18*20+11] = 7.109973e-02
	m[18*20+12] = 2.246759e-03
	m[18*20+13] = 7.074464e+00
	m[18*20+14] = 1.992133e-01
	m[18*20+15] = 8.104751e-01
	m[18*20+16] = 9.984255e-02
	m[18*20+17] = 6.121284e-01
	m[19*20+0] = 3.774477e+00
	m[19*20+1] = 1.366145e-01
	m[19*20+2] = 4.931206e-02
	m[19*20+3] = 4.076074e-01
	m[19*20+4] = 2.243512e-02
	m[19*20+5] = 9.047737e-03
	m[19*20+6] = 5.795409e-01
	m[19*20+7] = 4.228200e-01
	m[19*20+8] = 6.890244e+00
	m[19*20+9] = 7.926675e+00
	m[19*20+10] = 3.595310e+00
	m[19*20+11] = 3.493440e-02
	m[19*20+12] = 4.396720e+00
	m[19*20+13] = 1.643946e+00
	m[19*20+14] = 2.217442e-01
	m[19*20+15] = 7.477041e-02
	m[19*20+16] = 2.166054e-01
	m[19*20+17] = 9.663569e-02
	m[19*20+18] = 5.010635e-01

	for i = 0; i < naa; i++ {
		for j = 0; j < i; j++ {
			m[j*naa+i] = m[i*naa+j]
		}
	}

	pi[0] = 6.541704e-02
	pi[1] = 4.708366e-02
	pi[2] = 3.168984e-02
	pi[3] = 4.688141e-02
	pi[4] = 2.150693e-02
	pi[5] = 4.240711e-02
	pi[6] = 2.842211e-02
	pi[7] = 1.005278e-01
	pi[8] = 9.812606e-03
	pi[9] = 3.424424e-02
	pi[10] = 6.222565e-02
	pi[11] = 4.844488e-02
	pi[12] = 1.760370e-02
	pi[13] = 3.478555e-02
	pi[14] = 3.962469e-02
	pi[15] = 1.280566e-01
	pi[16] = 8.199314e-02
	pi[17] = 3.393045e-02
	pi[18] = 7.586119e-02
	pi[19] = 4.948141e-02

	dmat = mat.NewDense(naa, naa, m)
	return
}
