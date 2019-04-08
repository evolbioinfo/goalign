package protein

import (
	"github.com/gonum/matrix/mat64"
)

// Code taken from FastME

/*********************************************************/

/* Dayhoff's model data
 * Dayhoff, M.O., Schwartz, R.M., Orcutt, B.C. (1978)
 * "A model of evolutionary change in proteins."
 * Dayhoff, M.O.(ed.) Atlas of Protein Sequence Structur., Vol5, Suppl3.
 * National Biomedical Research Foundation, Washington DC, pp.345-352. */
func DayoffMats() (dmat *mat64.Dense, pi []float64) {
	var i, j, naa int
	naa = 20
	mat := make([]float64, naa*naa)
	pi = make([]float64, naa)

	mat[1*20+0] = 27.00
	mat[2*20+0] = 98.00
	mat[2*20+1] = 32.00
	mat[3*20+0] = 120.00
	mat[3*20+1] = 0.00
	mat[3*20+2] = 905.00
	mat[4*20+0] = 36.00
	mat[4*20+1] = 23.00
	mat[4*20+2] = 0.00
	mat[4*20+3] = 0.00
	mat[5*20+0] = 89.00
	mat[5*20+1] = 246.00
	mat[5*20+2] = 103.00
	mat[5*20+3] = 134.00
	mat[5*20+4] = 0.00
	mat[6*20+0] = 198.00
	mat[6*20+1] = 1.00
	mat[6*20+2] = 148.00
	mat[6*20+3] = 1153.00
	mat[6*20+4] = 0.00
	mat[6*20+5] = 716.00
	mat[7*20+0] = 240.00
	mat[7*20+1] = 9.00
	mat[7*20+2] = 139.00
	mat[7*20+3] = 125.00
	mat[7*20+4] = 11.00
	mat[7*20+5] = 28.00
	mat[7*20+6] = 81.00
	mat[8*20+0] = 23.00
	mat[8*20+1] = 240.00
	mat[8*20+2] = 535.00
	mat[8*20+3] = 86.00
	mat[8*20+4] = 28.00
	mat[8*20+5] = 606.00
	mat[8*20+6] = 43.00
	mat[8*20+7] = 10.00
	mat[9*20+0] = 65.00
	mat[9*20+1] = 64.00
	mat[9*20+2] = 77.00
	mat[9*20+3] = 24.00
	mat[9*20+4] = 44.00
	mat[9*20+5] = 18.00
	mat[9*20+6] = 61.00
	mat[9*20+7] = 0.00
	mat[9*20+8] = 7.00
	mat[10*20+0] = 41.00
	mat[10*20+1] = 15.00
	mat[10*20+2] = 34.00
	mat[10*20+3] = 0.00
	mat[10*20+4] = 0.00
	mat[10*20+5] = 73.00
	mat[10*20+6] = 11.00
	mat[10*20+7] = 7.00
	mat[10*20+8] = 44.00
	mat[10*20+9] = 257.00
	mat[11*20+0] = 26.00
	mat[11*20+1] = 464.00
	mat[11*20+2] = 318.00
	mat[11*20+3] = 71.00
	mat[11*20+4] = 0.00
	mat[11*20+5] = 153.00
	mat[11*20+6] = 83.00
	mat[11*20+7] = 27.00
	mat[11*20+8] = 26.00
	mat[11*20+9] = 46.00
	mat[11*20+10] = 18.00
	mat[12*20+0] = 72.00
	mat[12*20+1] = 90.00
	mat[12*20+2] = 1.00
	mat[12*20+3] = 0.00
	mat[12*20+4] = 0.00
	mat[12*20+5] = 114.00
	mat[12*20+6] = 30.00
	mat[12*20+7] = 17.00
	mat[12*20+8] = 0.00
	mat[12*20+9] = 336.00
	mat[12*20+10] = 527.00
	mat[12*20+11] = 243.00
	mat[13*20+0] = 18.00
	mat[13*20+1] = 14.00
	mat[13*20+2] = 14.00
	mat[13*20+3] = 0.00
	mat[13*20+4] = 0.00
	mat[13*20+5] = 0.00
	mat[13*20+6] = 0.00
	mat[13*20+7] = 15.00
	mat[13*20+8] = 48.00
	mat[13*20+9] = 196.00
	mat[13*20+10] = 157.00
	mat[13*20+11] = 0.00
	mat[13*20+12] = 92.00
	mat[14*20+0] = 250.00
	mat[14*20+1] = 103.00
	mat[14*20+2] = 42.00
	mat[14*20+3] = 13.00
	mat[14*20+4] = 19.00
	mat[14*20+5] = 153.00
	mat[14*20+6] = 51.00
	mat[14*20+7] = 34.00
	mat[14*20+8] = 94.00
	mat[14*20+9] = 12.00
	mat[14*20+10] = 32.00
	mat[14*20+11] = 33.00
	mat[14*20+12] = 17.00
	mat[14*20+13] = 11.00
	mat[15*20+0] = 409.00
	mat[15*20+1] = 154.00
	mat[15*20+2] = 495.00
	mat[15*20+3] = 95.00
	mat[15*20+4] = 161.00
	mat[15*20+5] = 56.00
	mat[15*20+6] = 79.00
	mat[15*20+7] = 234.00
	mat[15*20+8] = 35.00
	mat[15*20+9] = 24.00
	mat[15*20+10] = 17.00
	mat[15*20+11] = 96.00
	mat[15*20+12] = 62.00
	mat[15*20+13] = 46.00
	mat[15*20+14] = 245.00
	mat[16*20+0] = 371.00
	mat[16*20+1] = 26.00
	mat[16*20+2] = 229.00
	mat[16*20+3] = 66.00
	mat[16*20+4] = 16.00
	mat[16*20+5] = 53.00
	mat[16*20+6] = 34.00
	mat[16*20+7] = 30.00
	mat[16*20+8] = 22.00
	mat[16*20+9] = 192.00
	mat[16*20+10] = 33.00
	mat[16*20+11] = 136.00
	mat[16*20+12] = 104.00
	mat[16*20+13] = 13.00
	mat[16*20+14] = 78.00
	mat[16*20+15] = 550.00
	mat[17*20+0] = 0.00
	mat[17*20+1] = 201.00
	mat[17*20+2] = 23.00
	mat[17*20+3] = 0.00
	mat[17*20+4] = 0.00
	mat[17*20+5] = 0.00
	mat[17*20+6] = 0.00
	mat[17*20+7] = 0.00
	mat[17*20+8] = 27.00
	mat[17*20+9] = 0.00
	mat[17*20+10] = 46.00
	mat[17*20+11] = 0.00
	mat[17*20+12] = 0.00
	mat[17*20+13] = 76.00
	mat[17*20+14] = 0.00
	mat[17*20+15] = 75.00
	mat[17*20+16] = 0.00
	mat[18*20+0] = 24.00
	mat[18*20+1] = 8.00
	mat[18*20+2] = 95.00
	mat[18*20+3] = 0.00
	mat[18*20+4] = 96.00
	mat[18*20+5] = 0.00
	mat[18*20+6] = 22.00
	mat[18*20+7] = 0.00
	mat[18*20+8] = 127.00
	mat[18*20+9] = 37.00
	mat[18*20+10] = 28.00
	mat[18*20+11] = 13.00
	mat[18*20+12] = 0.00
	mat[18*20+13] = 698.00
	mat[18*20+14] = 0.00
	mat[18*20+15] = 34.00
	mat[18*20+16] = 42.00
	mat[18*20+17] = 61.00
	mat[19*20+0] = 208.00
	mat[19*20+1] = 24.00
	mat[19*20+2] = 15.00
	mat[19*20+3] = 18.00
	mat[19*20+4] = 49.00
	mat[19*20+5] = 35.00
	mat[19*20+6] = 37.00
	mat[19*20+7] = 54.00
	mat[19*20+8] = 44.00
	mat[19*20+9] = 889.00
	mat[19*20+10] = 175.00
	mat[19*20+11] = 10.00
	mat[19*20+12] = 258.00
	mat[19*20+13] = 12.00
	mat[19*20+14] = 48.00
	mat[19*20+15] = 30.00
	mat[19*20+16] = 157.00
	mat[19*20+17] = 0.00
	mat[19*20+18] = 28.00

	for i = 0; i < naa; i++ {
		for j = 0; j < i; j++ {
			mat[j*naa+i] = mat[i*naa+j]
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

	dmat = mat64.NewDense(naa, naa, mat)

	return
}

/*********************************************************/

/* JTT's model data
 * D.T.Jones, W.R.Taylor and J.M.Thornton
 * "The rapid generation of mutation data matrices from protein sequences"
 * CABIOS  vol.8 no.3 1992 pp275-282 */
func JTTMats() (dmat *mat64.Dense, pi []float64) {
	var i, j, naa int

	naa = 20

	mat := make([]float64, naa*naa)
	pi = make([]float64, naa)

	mat[1*20+0] = 58.00
	mat[2*20+0] = 54.00
	mat[2*20+1] = 45.00
	mat[3*20+0] = 81.00
	mat[3*20+1] = 16.00
	mat[3*20+2] = 528.00
	mat[4*20+0] = 56.00
	mat[4*20+1] = 113.00
	mat[4*20+2] = 34.00
	mat[4*20+3] = 10.00
	mat[5*20+0] = 57.00
	mat[5*20+1] = 310.00
	mat[5*20+2] = 86.00
	mat[5*20+3] = 49.00
	mat[5*20+4] = 9.00
	mat[6*20+0] = 105.00
	mat[6*20+1] = 29.00
	mat[6*20+2] = 58.00
	mat[6*20+3] = 767.00
	mat[6*20+4] = 5.00
	mat[6*20+5] = 323.00
	mat[7*20+0] = 179.00
	mat[7*20+1] = 137.00
	mat[7*20+2] = 81.00
	mat[7*20+3] = 130.00
	mat[7*20+4] = 59.00
	mat[7*20+5] = 26.00
	mat[7*20+6] = 119.00
	mat[8*20+0] = 27.00
	mat[8*20+1] = 328.00
	mat[8*20+2] = 391.00
	mat[8*20+3] = 112.00
	mat[8*20+4] = 69.00
	mat[8*20+5] = 597.00
	mat[8*20+6] = 26.00
	mat[8*20+7] = 23.00
	mat[9*20+0] = 36.00
	mat[9*20+1] = 22.00
	mat[9*20+2] = 47.00
	mat[9*20+3] = 11.00
	mat[9*20+4] = 17.00
	mat[9*20+5] = 9.00
	mat[9*20+6] = 12.00
	mat[9*20+7] = 6.00
	mat[9*20+8] = 16.00
	mat[10*20+0] = 30.00
	mat[10*20+1] = 38.00
	mat[10*20+2] = 12.00
	mat[10*20+3] = 7.00
	mat[10*20+4] = 23.00
	mat[10*20+5] = 72.00
	mat[10*20+6] = 9.00
	mat[10*20+7] = 6.00
	mat[10*20+8] = 56.00
	mat[10*20+9] = 229.00
	mat[11*20+0] = 35.00
	mat[11*20+1] = 646.00
	mat[11*20+2] = 263.00
	mat[11*20+3] = 26.00
	mat[11*20+4] = 7.00
	mat[11*20+5] = 292.00
	mat[11*20+6] = 181.00
	mat[11*20+7] = 27.00
	mat[11*20+8] = 45.00
	mat[11*20+9] = 21.00
	mat[11*20+10] = 14.00
	mat[12*20+0] = 54.00
	mat[12*20+1] = 44.00
	mat[12*20+2] = 30.00
	mat[12*20+3] = 15.00
	mat[12*20+4] = 31.00
	mat[12*20+5] = 43.00
	mat[12*20+6] = 18.00
	mat[12*20+7] = 14.00
	mat[12*20+8] = 33.00
	mat[12*20+9] = 479.00
	mat[12*20+10] = 388.00
	mat[12*20+11] = 65.00
	mat[13*20+0] = 15.00
	mat[13*20+1] = 5.00
	mat[13*20+2] = 10.00
	mat[13*20+3] = 4.00
	mat[13*20+4] = 78.00
	mat[13*20+5] = 4.00
	mat[13*20+6] = 5.00
	mat[13*20+7] = 5.00
	mat[13*20+8] = 40.00
	mat[13*20+9] = 89.00
	mat[13*20+10] = 248.00
	mat[13*20+11] = 4.00
	mat[13*20+12] = 43.00
	mat[14*20+0] = 194.00
	mat[14*20+1] = 74.00
	mat[14*20+2] = 15.00
	mat[14*20+3] = 15.00
	mat[14*20+4] = 14.00
	mat[14*20+5] = 164.00
	mat[14*20+6] = 18.00
	mat[14*20+7] = 24.00
	mat[14*20+8] = 115.00
	mat[14*20+9] = 10.00
	mat[14*20+10] = 102.00
	mat[14*20+11] = 21.00
	mat[14*20+12] = 16.00
	mat[14*20+13] = 17.00
	mat[15*20+0] = 378.00
	mat[15*20+1] = 101.00
	mat[15*20+2] = 503.00
	mat[15*20+3] = 59.00
	mat[15*20+4] = 223.00
	mat[15*20+5] = 53.00
	mat[15*20+6] = 30.00
	mat[15*20+7] = 201.00
	mat[15*20+8] = 73.00
	mat[15*20+9] = 40.00
	mat[15*20+10] = 59.00
	mat[15*20+11] = 47.00
	mat[15*20+12] = 29.00
	mat[15*20+13] = 92.00
	mat[15*20+14] = 285.00
	mat[16*20+0] = 475.00
	mat[16*20+1] = 64.00
	mat[16*20+2] = 232.00
	mat[16*20+3] = 38.00
	mat[16*20+4] = 42.00
	mat[16*20+5] = 51.00
	mat[16*20+6] = 32.00
	mat[16*20+7] = 33.00
	mat[16*20+8] = 46.00
	mat[16*20+9] = 245.00
	mat[16*20+10] = 25.00
	mat[16*20+11] = 103.00
	mat[16*20+12] = 226.00
	mat[16*20+13] = 12.00
	mat[16*20+14] = 118.00
	mat[16*20+15] = 477.00
	mat[17*20+0] = 9.00
	mat[17*20+1] = 126.00
	mat[17*20+2] = 8.00
	mat[17*20+3] = 4.00
	mat[17*20+4] = 115.00
	mat[17*20+5] = 18.00
	mat[17*20+6] = 10.00
	mat[17*20+7] = 55.00
	mat[17*20+8] = 8.00
	mat[17*20+9] = 9.00
	mat[17*20+10] = 52.00
	mat[17*20+11] = 10.00
	mat[17*20+12] = 24.00
	mat[17*20+13] = 53.00
	mat[17*20+14] = 6.00
	mat[17*20+15] = 35.00
	mat[17*20+16] = 12.00
	mat[18*20+0] = 11.00
	mat[18*20+1] = 20.00
	mat[18*20+2] = 70.00
	mat[18*20+3] = 46.00
	mat[18*20+4] = 209.00
	mat[18*20+5] = 24.00
	mat[18*20+6] = 7.00
	mat[18*20+7] = 8.00
	mat[18*20+8] = 573.00
	mat[18*20+9] = 32.00
	mat[18*20+10] = 24.00
	mat[18*20+11] = 8.00
	mat[18*20+12] = 18.00
	mat[18*20+13] = 536.00
	mat[18*20+14] = 10.00
	mat[18*20+15] = 63.00
	mat[18*20+16] = 21.00
	mat[18*20+17] = 71.00
	mat[19*20+0] = 298.00
	mat[19*20+1] = 17.00
	mat[19*20+2] = 16.00
	mat[19*20+3] = 31.00
	mat[19*20+4] = 62.00
	mat[19*20+5] = 20.00
	mat[19*20+6] = 45.00
	mat[19*20+7] = 47.00
	mat[19*20+8] = 11.00
	mat[19*20+9] = 961.00
	mat[19*20+10] = 180.00
	mat[19*20+11] = 14.00
	mat[19*20+12] = 323.00
	mat[19*20+13] = 62.00
	mat[19*20+14] = 23.00
	mat[19*20+15] = 38.00
	mat[19*20+16] = 112.00
	mat[19*20+17] = 25.00
	mat[19*20+18] = 16.00

	for i = 0; i < naa; i++ {
		for j = 0; j < i; j++ {
			mat[j*naa+i] = mat[i*naa+j]
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

	dmat = mat64.NewDense(naa, naa, mat)
	return
}

/*********************************************************/

func MtREVMats() (dmat *mat64.Dense, pi []float64) {
	var i, j, naa int

	naa = 20

	mat := make([]float64, naa*naa)
	pi = make([]float64, naa)

	mat[1*20+0] = 23.18
	mat[2*20+0] = 26.95
	mat[2*20+1] = 13.24
	mat[3*20+0] = 17.67
	mat[3*20+1] = 1.90
	mat[3*20+2] = 794.38
	mat[4*20+0] = 59.93
	mat[4*20+1] = 103.33
	mat[4*20+2] = 58.94
	mat[4*20+3] = 1.90
	mat[5*20+0] = 1.90
	mat[5*20+1] = 220.99
	mat[5*20+2] = 173.56
	mat[5*20+3] = 55.28
	mat[5*20+4] = 75.24
	mat[6*20+0] = 9.77
	mat[6*20+1] = 1.90
	mat[6*20+2] = 63.05
	mat[6*20+3] = 583.55
	mat[6*20+4] = 1.90
	mat[6*20+5] = 313.56
	mat[7*20+0] = 120.71
	mat[7*20+1] = 23.03
	mat[7*20+2] = 53.30
	mat[7*20+3] = 56.77
	mat[7*20+4] = 30.71
	mat[7*20+5] = 6.75
	mat[7*20+6] = 28.28
	mat[8*20+0] = 13.90
	mat[8*20+1] = 165.23
	mat[8*20+2] = 496.13
	mat[8*20+3] = 113.99
	mat[8*20+4] = 141.49
	mat[8*20+5] = 582.40
	mat[8*20+6] = 49.12
	mat[8*20+7] = 1.90
	mat[9*20+0] = 96.49
	mat[9*20+1] = 1.90
	mat[9*20+2] = 27.10
	mat[9*20+3] = 4.34
	mat[9*20+4] = 62.73
	mat[9*20+5] = 8.34
	mat[9*20+6] = 3.31
	mat[9*20+7] = 5.98
	mat[9*20+8] = 12.26
	mat[10*20+0] = 25.46
	mat[10*20+1] = 15.58
	mat[10*20+2] = 15.16
	mat[10*20+3] = 1.90
	mat[10*20+4] = 25.65
	mat[10*20+5] = 39.70
	mat[10*20+6] = 1.90
	mat[10*20+7] = 2.41
	mat[10*20+8] = 11.49
	mat[10*20+9] = 329.09
	mat[11*20+0] = 8.36
	mat[11*20+1] = 141.40
	mat[11*20+2] = 608.70
	mat[11*20+3] = 2.31
	mat[11*20+4] = 1.90
	mat[11*20+5] = 465.58
	mat[11*20+6] = 313.86
	mat[11*20+7] = 22.73
	mat[11*20+8] = 127.67
	mat[11*20+9] = 19.57
	mat[11*20+10] = 14.88
	mat[12*20+0] = 141.88
	mat[12*20+1] = 1.90
	mat[12*20+2] = 65.41
	mat[12*20+3] = 1.90
	mat[12*20+4] = 6.18
	mat[12*20+5] = 47.37
	mat[12*20+6] = 1.90
	mat[12*20+7] = 1.90
	mat[12*20+8] = 11.97
	mat[12*20+9] = 517.98
	mat[12*20+10] = 537.53
	mat[12*20+11] = 91.37
	mat[13*20+0] = 6.37
	mat[13*20+1] = 4.69
	mat[13*20+2] = 15.20
	mat[13*20+3] = 4.98
	mat[13*20+4] = 70.80
	mat[13*20+5] = 19.11
	mat[13*20+6] = 2.67
	mat[13*20+7] = 1.90
	mat[13*20+8] = 48.16
	mat[13*20+9] = 84.67
	mat[13*20+10] = 216.06
	mat[13*20+11] = 6.44
	mat[13*20+12] = 90.82
	mat[14*20+0] = 54.31
	mat[14*20+1] = 23.64
	mat[14*20+2] = 73.31
	mat[14*20+3] = 13.43
	mat[14*20+4] = 31.26
	mat[14*20+5] = 137.29
	mat[14*20+6] = 12.83
	mat[14*20+7] = 1.90
	mat[14*20+8] = 60.97
	mat[14*20+9] = 20.63
	mat[14*20+10] = 40.10
	mat[14*20+11] = 50.10
	mat[14*20+12] = 18.84
	mat[14*20+13] = 17.31
	mat[15*20+0] = 387.86
	mat[15*20+1] = 6.04
	mat[15*20+2] = 494.39
	mat[15*20+3] = 69.02
	mat[15*20+4] = 277.05
	mat[15*20+5] = 54.11
	mat[15*20+6] = 54.71
	mat[15*20+7] = 125.93
	mat[15*20+8] = 77.46
	mat[15*20+9] = 47.70
	mat[15*20+10] = 73.61
	mat[15*20+11] = 105.79
	mat[15*20+12] = 111.16
	mat[15*20+13] = 64.29
	mat[15*20+14] = 169.90
	mat[16*20+0] = 480.72
	mat[16*20+1] = 2.08
	mat[16*20+2] = 238.46
	mat[16*20+3] = 28.01
	mat[16*20+4] = 179.97
	mat[16*20+5] = 94.93
	mat[16*20+6] = 14.82
	mat[16*20+7] = 11.17
	mat[16*20+8] = 44.78
	mat[16*20+9] = 368.43
	mat[16*20+10] = 126.40
	mat[16*20+11] = 136.33
	mat[16*20+12] = 528.17
	mat[16*20+13] = 33.85
	mat[16*20+14] = 128.22
	mat[16*20+15] = 597.21
	mat[17*20+0] = 1.90
	mat[17*20+1] = 21.95
	mat[17*20+2] = 10.68
	mat[17*20+3] = 19.86
	mat[17*20+4] = 33.60
	mat[17*20+5] = 1.90
	mat[17*20+6] = 1.90
	mat[17*20+7] = 10.92
	mat[17*20+8] = 7.08
	mat[17*20+9] = 1.90
	mat[17*20+10] = 32.44
	mat[17*20+11] = 24.00
	mat[17*20+12] = 21.71
	mat[17*20+13] = 7.84
	mat[17*20+14] = 4.21
	mat[17*20+15] = 38.58
	mat[17*20+16] = 9.99
	mat[18*20+0] = 6.48
	mat[18*20+1] = 1.90
	mat[18*20+2] = 191.36
	mat[18*20+3] = 21.21
	mat[18*20+4] = 254.77
	mat[18*20+5] = 38.82
	mat[18*20+6] = 13.12
	mat[18*20+7] = 3.21
	mat[18*20+8] = 670.14
	mat[18*20+9] = 25.01
	mat[18*20+10] = 44.15
	mat[18*20+11] = 51.17
	mat[18*20+12] = 39.96
	mat[18*20+13] = 465.58
	mat[18*20+14] = 16.21
	mat[18*20+15] = 64.92
	mat[18*20+16] = 38.73
	mat[18*20+17] = 26.25
	mat[19*20+0] = 195.06
	mat[19*20+1] = 7.64
	mat[19*20+2] = 1.90
	mat[19*20+3] = 1.90
	mat[19*20+4] = 1.90
	mat[19*20+5] = 19.00
	mat[19*20+6] = 21.14
	mat[19*20+7] = 2.53
	mat[19*20+8] = 1.90
	mat[19*20+9] = 1222.94
	mat[19*20+10] = 91.67
	mat[19*20+11] = 1.90
	mat[19*20+12] = 387.54
	mat[19*20+13] = 6.35
	mat[19*20+14] = 8.23
	mat[19*20+15] = 1.90
	mat[19*20+16] = 204.54
	mat[19*20+17] = 5.37
	mat[19*20+18] = 1.90

	for i = 0; i < naa; i++ {
		for j = 0; j < i; j++ {
			mat[j*naa+i] = mat[i*naa+j]
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

	dmat = mat64.NewDense(naa, naa, mat)
	return
}

/*********************************************************/

/* LG model
 * Si Quang LE & Olivier Gascuel
 * "An improved general amino-acid replacement matrix"
 * Mol Biol Evol. 2008 Jul;25(7):1307-20. */
func LGMats() (dmat *mat64.Dense, pi []float64) {
	var i, j, naa int

	naa = 20

	mat := make([]float64, naa*naa)
	pi = make([]float64, naa)

	mat[1*20+0] = 0.449682
	mat[2*20+0] = 0.267582
	mat[2*20+1] = 0.827348
	mat[3*20+0] = 0.401081
	mat[3*20+1] = 0.132811
	mat[3*20+2] = 5.921004
	mat[4*20+0] = 2.312843
	mat[4*20+1] = 0.552587
	mat[4*20+2] = 0.522133
	mat[4*20+3] = 0.056428
	mat[5*20+0] = 0.944706
	mat[5*20+1] = 3.109412
	mat[5*20+2] = 1.877436
	mat[5*20+3] = 0.498202
	mat[5*20+4] = 0.080602
	mat[6*20+0] = 1.164358
	mat[6*20+1] = 0.442407
	mat[6*20+2] = 0.599223
	mat[6*20+3] = 6.374225
	mat[6*20+4] = 0.001330
	mat[6*20+5] = 4.799804
	mat[7*20+0] = 2.101845
	mat[7*20+1] = 0.443980
	mat[7*20+2] = 1.566189
	mat[7*20+3] = 0.922928
	mat[7*20+4] = 0.529114
	mat[7*20+5] = 0.279365
	mat[7*20+6] = 0.407773
	mat[8*20+0] = 0.341479
	mat[8*20+1] = 2.657648
	mat[8*20+2] = 4.889564
	mat[8*20+3] = 0.982202
	mat[8*20+4] = 0.593147
	mat[8*20+5] = 5.177996
	mat[8*20+6] = 0.458209
	mat[8*20+7] = 0.304320
	mat[9*20+0] = 0.122945
	mat[9*20+1] = 0.134451
	mat[9*20+2] = 0.216069
	mat[9*20+3] = 0.010922
	mat[9*20+4] = 0.262931
	mat[9*20+5] = 0.073719
	mat[9*20+6] = 0.056153
	mat[9*20+7] = 0.008454
	mat[9*20+8] = 0.106232
	mat[10*20+0] = 0.391826
	mat[10*20+1] = 0.330360
	mat[10*20+2] = 0.075149
	mat[10*20+3] = 0.017176
	mat[10*20+4] = 0.541544
	mat[10*20+5] = 0.613290
	mat[10*20+6] = 0.086633
	mat[10*20+7] = 0.047556
	mat[10*20+8] = 0.363554
	mat[10*20+9] = 3.801506
	mat[11*20+0] = 0.556137
	mat[11*20+1] = 7.114371
	mat[11*20+2] = 2.463341
	mat[11*20+3] = 0.278545
	mat[11*20+4] = 0.003892
	mat[11*20+5] = 3.466773
	mat[11*20+6] = 2.168935
	mat[11*20+7] = 0.313114
	mat[11*20+8] = 0.682564
	mat[11*20+9] = 0.173179
	mat[11*20+10] = 0.145273
	mat[12*20+0] = 1.050301
	mat[12*20+1] = 0.477124
	mat[12*20+2] = 0.370061
	mat[12*20+3] = 0.022762
	mat[12*20+4] = 0.773189
	mat[12*20+5] = 1.656669
	mat[12*20+6] = 0.183748
	mat[12*20+7] = 0.137976
	mat[12*20+8] = 0.395265
	mat[12*20+9] = 3.849020
	mat[12*20+10] = 5.836269
	mat[12*20+11] = 0.672252
	mat[13*20+0] = 0.237746
	mat[13*20+1] = 0.055544
	mat[13*20+2] = 0.090929
	mat[13*20+3] = 0.017714
	mat[13*20+4] = 0.950511
	mat[13*20+5] = 0.033627
	mat[13*20+6] = 0.024362
	mat[13*20+7] = 0.080743
	mat[13*20+8] = 0.616582
	mat[13*20+9] = 1.020659
	mat[13*20+10] = 2.426267
	mat[13*20+11] = 0.026721
	mat[13*20+12] = 1.626175
	mat[14*20+0] = 1.232907
	mat[14*20+1] = 0.404818
	mat[14*20+2] = 0.190630
	mat[14*20+3] = 0.449817
	mat[14*20+4] = 0.076565
	mat[14*20+5] = 0.698390
	mat[14*20+6] = 0.523437
	mat[14*20+7] = 0.226307
	mat[14*20+8] = 0.545492
	mat[14*20+9] = 0.086269
	mat[14*20+10] = 0.265077
	mat[14*20+11] = 0.445474
	mat[14*20+12] = 0.096861
	mat[14*20+13] = 0.104849
	mat[15*20+0] = 4.655234
	mat[15*20+1] = 0.897892
	mat[15*20+2] = 4.299421
	mat[15*20+3] = 1.268215
	mat[15*20+4] = 2.605967
	mat[15*20+5] = 1.205796
	mat[15*20+6] = 0.667092
	mat[15*20+7] = 1.784779
	mat[15*20+8] = 0.947402
	mat[15*20+9] = 0.063251
	mat[15*20+10] = 0.184361
	mat[15*20+11] = 0.755746
	mat[15*20+12] = 0.319101
	mat[15*20+13] = 0.355654
	mat[15*20+14] = 1.424806
	mat[16*20+0] = 1.986433
	mat[16*20+1] = 0.579784
	mat[16*20+2] = 2.061491
	mat[16*20+3] = 0.405969
	mat[16*20+4] = 0.993542
	mat[16*20+5] = 1.027335
	mat[16*20+6] = 0.659097
	mat[16*20+7] = 0.114336
	mat[16*20+8] = 0.526423
	mat[16*20+9] = 0.992803
	mat[16*20+10] = 0.286481
	mat[16*20+11] = 1.152184
	mat[16*20+12] = 1.866946
	mat[16*20+13] = 0.145526
	mat[16*20+14] = 0.592443
	mat[16*20+15] = 6.266071
	mat[17*20+0] = 0.179433
	mat[17*20+1] = 0.701255
	mat[17*20+2] = 0.054722
	mat[17*20+3] = 0.046559
	mat[17*20+4] = 0.659458
	mat[17*20+5] = 0.249044
	mat[17*20+6] = 0.099542
	mat[17*20+7] = 0.292882
	mat[17*20+8] = 0.559689
	mat[17*20+9] = 0.121839
	mat[17*20+10] = 0.649934
	mat[17*20+11] = 0.047995
	mat[17*20+12] = 0.660667
	mat[17*20+13] = 2.425821
	mat[17*20+14] = 0.118287
	mat[17*20+15] = 0.267487
	mat[17*20+16] = 0.144967
	mat[18*20+0] = 0.223517
	mat[18*20+1] = 0.342216
	mat[18*20+2] = 0.658002
	mat[18*20+3] = 0.147235
	mat[18*20+4] = 1.095311
	mat[18*20+5] = 0.244886
	mat[18*20+6] = 0.140547
	mat[18*20+7] = 0.056885
	mat[18*20+8] = 5.446234
	mat[18*20+9] = 0.238891
	mat[18*20+10] = 0.292232
	mat[18*20+11] = 0.138336
	mat[18*20+12] = 0.436403
	mat[18*20+13] = 7.598781
	mat[18*20+14] = 0.109774
	mat[18*20+15] = 0.407468
	mat[18*20+16] = 0.236493
	mat[18*20+17] = 3.344523
	mat[19*20+0] = 2.368823
	mat[19*20+1] = 0.173721
	mat[19*20+2] = 0.088856
	mat[19*20+3] = 0.038720
	mat[19*20+4] = 1.745884
	mat[19*20+5] = 0.204644
	mat[19*20+6] = 0.278624
	mat[19*20+7] = 0.075577
	mat[19*20+8] = 0.108961
	mat[19*20+9] = 9.416771
	mat[19*20+10] = 1.519645
	mat[19*20+11] = 0.184432
	mat[19*20+12] = 1.595049
	mat[19*20+13] = 0.578417
	mat[19*20+14] = 0.302548
	mat[19*20+15] = 0.062285
	mat[19*20+16] = 1.947321
	mat[19*20+17] = 0.201078
	mat[19*20+18] = 0.235819

	for i = 0; i < naa; i++ {
		for j = 0; j < i; j++ {
			mat[j*naa+i] = mat[i*naa+j]
		}
	}
	pi[0] = 0.079611
	pi[1] = 0.053191
	pi[2] = 0.039948
	pi[3] = 0.050634
	pi[4] = 0.013590
	pi[5] = 0.038611
	pi[6] = 0.066539
	pi[7] = 0.059913
	pi[8] = 0.021738
	pi[9] = 0.063589
	pi[10] = 0.105134
	pi[11] = 0.061845
	pi[12] = 0.022990
	pi[13] = 0.044365
	pi[14] = 0.044909
	pi[15] = 0.059477
	pi[16] = 0.054114
	pi[17] = 0.012588
	pi[18] = 0.035709
	pi[19] = 0.071505

	dmat = mat64.NewDense(naa, naa, mat)
	return
}

/*********************************************************/

/* WAG's model data
 * Simon Whelan and Nick Goldman
 * "A general empirical model of protein evolution derived from multiple
 * protein families using a maximum-likelihood approach"
 * MBE (2001) 18:691-699 */
func WAGMats() (dmat *mat64.Dense, pi []float64) {
	var i, j, naa int

	naa = 20

	mat := make([]float64, naa*naa)
	pi = make([]float64, naa)

	mat[1*20+0] = 55.15710
	mat[2*20+0] = 50.98480
	mat[2*20+1] = 63.53460
	mat[3*20+0] = 73.89980
	mat[3*20+1] = 14.73040
	mat[3*20+2] = 542.94200
	mat[4*20+0] = 102.70400
	mat[4*20+1] = 52.81910
	mat[4*20+2] = 26.52560
	mat[4*20+3] = 3.02949
	mat[5*20+0] = 90.85980
	mat[5*20+1] = 303.55000
	mat[5*20+2] = 154.36400
	mat[5*20+3] = 61.67830
	mat[5*20+4] = 9.88179
	mat[6*20+0] = 158.28500
	mat[6*20+1] = 43.91570
	mat[6*20+2] = 94.71980
	mat[6*20+3] = 617.41600
	mat[6*20+4] = 2.13520
	mat[6*20+5] = 546.94700
	mat[7*20+0] = 141.67200
	mat[7*20+1] = 58.46650
	mat[7*20+2] = 112.55600
	mat[7*20+3] = 86.55840
	mat[7*20+4] = 30.66740
	mat[7*20+5] = 33.00520
	mat[7*20+6] = 56.77170
	mat[8*20+0] = 31.69540
	mat[8*20+1] = 213.71500
	mat[8*20+2] = 395.62900
	mat[8*20+3] = 93.06760
	mat[8*20+4] = 24.89720
	mat[8*20+5] = 429.41100
	mat[8*20+6] = 57.00250
	mat[8*20+7] = 24.94100
	mat[9*20+0] = 19.33350
	mat[9*20+1] = 18.69790
	mat[9*20+2] = 55.42360
	mat[9*20+3] = 3.94370
	mat[9*20+4] = 17.01350
	mat[9*20+5] = 11.39170
	mat[9*20+6] = 12.73950
	mat[9*20+7] = 3.04501
	mat[9*20+8] = 13.81900
	mat[10*20+0] = 39.79150
	mat[10*20+1] = 49.76710
	mat[10*20+2] = 13.15280
	mat[10*20+3] = 8.48047
	mat[10*20+4] = 38.42870
	mat[10*20+5] = 86.94890
	mat[10*20+6] = 15.42630
	mat[10*20+7] = 6.13037
	mat[10*20+8] = 49.94620
	mat[10*20+9] = 317.09700
	mat[11*20+0] = 90.62650
	mat[11*20+1] = 535.14200
	mat[11*20+2] = 301.20100
	mat[11*20+3] = 47.98550
	mat[11*20+4] = 7.40339
	mat[11*20+5] = 389.49000
	mat[11*20+6] = 258.44300
	mat[11*20+7] = 37.35580
	mat[11*20+8] = 89.04320
	mat[11*20+9] = 32.38320
	mat[11*20+10] = 25.75550
	mat[12*20+0] = 89.34960
	mat[12*20+1] = 68.31620
	mat[12*20+2] = 19.82210
	mat[12*20+3] = 10.37540
	mat[12*20+4] = 39.04820
	mat[12*20+5] = 154.52600
	mat[12*20+6] = 31.51240
	mat[12*20+7] = 17.41000
	mat[12*20+8] = 40.41410
	mat[12*20+9] = 425.74600
	mat[12*20+10] = 485.40200
	mat[12*20+11] = 93.42760
	mat[13*20+0] = 21.04940
	mat[13*20+1] = 10.27110
	mat[13*20+2] = 9.61621
	mat[13*20+3] = 4.67304
	mat[13*20+4] = 39.80200
	mat[13*20+5] = 9.99208
	mat[13*20+6] = 8.11339
	mat[13*20+7] = 4.99310
	mat[13*20+8] = 67.93710
	mat[13*20+9] = 105.94700
	mat[13*20+10] = 211.51700
	mat[13*20+11] = 8.88360
	mat[13*20+12] = 119.06300
	mat[14*20+0] = 143.85500
	mat[14*20+1] = 67.94890
	mat[14*20+2] = 19.50810
	mat[14*20+3] = 42.39840
	mat[14*20+4] = 10.94040
	mat[14*20+5] = 93.33720
	mat[14*20+6] = 68.23550
	mat[14*20+7] = 24.35700
	mat[14*20+8] = 69.61980
	mat[14*20+9] = 9.99288
	mat[14*20+10] = 41.58440
	mat[14*20+11] = 55.68960
	mat[14*20+12] = 17.13290
	mat[14*20+13] = 16.14440
	mat[15*20+0] = 337.07900
	mat[15*20+1] = 122.41900
	mat[15*20+2] = 397.42300
	mat[15*20+3] = 107.17600
	mat[15*20+4] = 140.76600
	mat[15*20+5] = 102.88700
	mat[15*20+6] = 70.49390
	mat[15*20+7] = 134.18200
	mat[15*20+8] = 74.01690
	mat[15*20+9] = 31.94400
	mat[15*20+10] = 34.47390
	mat[15*20+11] = 96.71300
	mat[15*20+12] = 49.39050
	mat[15*20+13] = 54.59310
	mat[15*20+14] = 161.32800
	mat[16*20+0] = 212.11100
	mat[16*20+1] = 55.44130
	mat[16*20+2] = 203.00600
	mat[16*20+3] = 37.48660
	mat[16*20+4] = 51.29840
	mat[16*20+5] = 85.79280
	mat[16*20+6] = 82.27650
	mat[16*20+7] = 22.58330
	mat[16*20+8] = 47.33070
	mat[16*20+9] = 145.81600
	mat[16*20+10] = 32.66220
	mat[16*20+11] = 138.69800
	mat[16*20+12] = 151.61200
	mat[16*20+13] = 17.19030
	mat[16*20+14] = 79.53840
	mat[16*20+15] = 437.80200
	mat[17*20+0] = 11.31330
	mat[17*20+1] = 116.39200
	mat[17*20+2] = 7.19167
	mat[17*20+3] = 12.97670
	mat[17*20+4] = 71.70700
	mat[17*20+5] = 21.57370
	mat[17*20+6] = 15.65570
	mat[17*20+7] = 33.69830
	mat[17*20+8] = 26.25690
	mat[17*20+9] = 21.24830
	mat[17*20+10] = 66.53090
	mat[17*20+11] = 13.75050
	mat[17*20+12] = 51.57060
	mat[17*20+13] = 152.96400
	mat[17*20+14] = 13.94050
	mat[17*20+15] = 52.37420
	mat[17*20+16] = 11.08640
	mat[18*20+0] = 24.07350
	mat[18*20+1] = 38.15330
	mat[18*20+2] = 108.60000
	mat[18*20+3] = 32.57110
	mat[18*20+4] = 54.38330
	mat[18*20+5] = 22.77100
	mat[18*20+6] = 19.63030
	mat[18*20+7] = 10.36040
	mat[18*20+8] = 387.34400
	mat[18*20+9] = 42.01700
	mat[18*20+10] = 39.86180
	mat[18*20+11] = 13.32640
	mat[18*20+12] = 42.84370
	mat[18*20+13] = 645.42800
	mat[18*20+14] = 21.60460
	mat[18*20+15] = 78.69930
	mat[18*20+16] = 29.11480
	mat[18*20+17] = 248.53900
	mat[19*20+0] = 200.60100
	mat[19*20+1] = 25.18490
	mat[19*20+2] = 19.62460
	mat[19*20+3] = 15.23350
	mat[19*20+4] = 100.21400
	mat[19*20+5] = 30.12810
	mat[19*20+6] = 58.87310
	mat[19*20+7] = 18.72470
	mat[19*20+8] = 11.83580
	mat[19*20+9] = 782.13000
	mat[19*20+10] = 180.03400
	mat[19*20+11] = 30.54340
	mat[19*20+12] = 205.84500
	mat[19*20+13] = 64.98920
	mat[19*20+14] = 31.48870
	mat[19*20+15] = 23.27390
	mat[19*20+16] = 138.82300
	mat[19*20+17] = 36.53690
	mat[19*20+18] = 31.47300

	for i = 0; i < naa; i++ {
		for j = 0; j < i; j++ {
			mat[j*naa+i] = mat[i*naa+j]
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

	dmat = mat64.NewDense(naa, naa, mat)
	return
}
