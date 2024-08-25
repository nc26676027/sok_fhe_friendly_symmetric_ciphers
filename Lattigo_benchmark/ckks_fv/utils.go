package ckks_fv

import (
	"encoding/binary"
	"io"
	"math"
	"math/big"
	"math/bits"

	"github.com/ldsec/lattigo/v2/ring"
)

func htobe64(value uint64) []byte {
	buf := make([]byte, 8)
	binary.BigEndian.PutUint64(buf, value)
	return buf
}

func calculateRow(prevRow, firstRow []uint64, q uint64) (out []uint64){
	
	PASTA_T := len(prevRow)
	out = make([]uint64, PASTA_T)
	for i := 0; i < PASTA_T; i++{
		tmp := new(big.Int).Mul(new(big.Int).SetUint64(firstRow[i]), new(big.Int).SetUint64(prevRow[PASTA_T-1]))
		tmp.Mod(tmp, new(big.Int).SetUint64(q) )
		if i > 0 {
			tmp.Add(tmp, new(big.Int).SetUint64(prevRow[i-1]))
			tmp.Mod(tmp, new(big.Int).SetUint64(q) )
		}
		out[i] = tmp.Uint64()
	}
	return out
}

func GetRandomVector( vecSize int, rand io.Reader, q uint64 ) (vec []uint64){
	vec = make([]uint64, vecSize)
	for i := 0; i < vecSize; i++ {
		vec[i] = SampleZqx(rand, q)
	}
	return vec
}

//----------------------------------------------------------------
// (0  1  0  ... 0 ) ^ t
// (0  0  1  ... 0 )
// (.  .  .  ... . )
// (.  .  .  ... . )
// (.  .  .  ... . )
// (0  0  0 ...  1 )
// (r1 r2 r3 ... rt)
// invertible by design
func GetRandomMatrixPasta( matSize int, rand io.Reader, q uint64) (mat [][]uint64){
	mat = make([][]uint64, matSize)
	mat[0] = GetRandomVector(matSize, rand, q)
	for i := 1; i < matSize; i++ {
		mat[i] = calculateRow(mat[i-1], mat[0], q)
	}  
	return mat
}

func GetRandomMatrixMasta( matSize int, rand io.Reader, q uint64) (mat [][]uint64){
	MastaAlpha := 3

	mat = make([][]uint64, matSize)
	for i := 0; i < len(mat); i++ {
		mat[i] = make([]uint64, matSize)
	}

	fieldElement := GetRandomVector(matSize, rand, q)
	for col := 0; col < matSize; col++ {
		for row := 0; row < matSize; row++ {
			ind := (row - col)
			if ind < 0 {
				tmp := new(big.Int).Mul( new(big.Int).SetUint64(fieldElement[ind + matSize]), new(big.Int).SetUint64(uint64(MastaAlpha)) )
				tmp.Mod(tmp, new(big.Int).SetUint64(q) )
				mat[row][col] = tmp.Uint64()
			} else {
				mat[row][col] = fieldElement[ind]
			}
		}
	}
	return mat
}

// Returns uniform random value in (0,q) by rejection sampling
func SampleZqx(rand io.Reader, q uint64) (res uint64) {
	bitLen := bits.Len64(q - 2)
	byteLen := (bitLen + 7) / 8
	b := bitLen % 8
	if b == 0 {
		b = 8
	}

	bytes := make([]byte, byteLen)
	for {
		_, err := io.ReadFull(rand, bytes)
		if err != nil {
			panic(err)
		}
		bytes[byteLen-1] &= uint8((1 << b) - 1)

		res = 0
		for i := 0; i < byteLen; i++ {
			res += uint64(bytes[i]) << (8 * i)
		}

		if res < q {
			return
		}
	}
}

// StandardDeviation computes the scaled standard deviation of the input vector.
func StandardDeviation(vec []float64, scale float64) (std float64) {
	// We assume that the error is centered around zero
	var err, tmp, mean, n float64

	n = float64(len(vec))

	for _, c := range vec {
		mean += c
	}

	mean /= n

	for _, c := range vec {
		tmp = c - mean
		err += tmp * tmp
	}

	return math.Sqrt(err/n) * scale
}

func scaleUpExact(value float64, n float64, q uint64) (res uint64) {

	var isNegative bool
	var xFlo *big.Float
	var xInt *big.Int

	isNegative = false
	if value < 0 {
		isNegative = true
		xFlo = big.NewFloat(-n * value)
	} else {
		xFlo = big.NewFloat(n * value)
	}

	xFlo.Add(xFlo, big.NewFloat(0.5))

	xInt = new(big.Int)
	xFlo.Int(xInt)
	xInt.Mod(xInt, ring.NewUint(q))

	res = xInt.Uint64()

	if isNegative {
		res = q - res
	}

	return
}

func scaleUpVecExact(values []float64, n float64, moduli []uint64, coeffs [][]uint64) {

	var isNegative bool
	var xFlo *big.Float
	var xInt *big.Int
	tmp := new(big.Int)

	for i := range values {

		if n*math.Abs(values[i]) > 1.8446744073709552e+19 {

			isNegative = false
			if values[i] < 0 {
				isNegative = true
				xFlo = big.NewFloat(-n * values[i])
			} else {
				xFlo = big.NewFloat(n * values[i])
			}

			xFlo.Add(xFlo, big.NewFloat(0.5))

			xInt = new(big.Int)
			xFlo.Int(xInt)

			for j := range moduli {
				tmp.Mod(xInt, ring.NewUint(moduli[j]))
				if isNegative {
					coeffs[j][i] = moduli[j] - tmp.Uint64()
				} else {
					coeffs[j][i] = tmp.Uint64()
				}
			}
		} else {

			if values[i] < 0 {
				for j := range moduli {
					coeffs[j][i] = moduli[j] - (uint64(-n*values[i]+0.5) % moduli[j])
				}
			} else {
				for j := range moduli {
					coeffs[j][i] = uint64(n*values[i]+0.5) % moduli[j]
				}
			}
		}
	}
}

func scaleUpVecExactBigFloat(values []*big.Float, scale float64, moduli []uint64, coeffs [][]uint64) {

	prec := int(values[0].Prec())

	xFlo := ring.NewFloat(0, prec)
	xInt := new(big.Int)
	tmp := new(big.Int)

	zero := ring.NewFloat(0, prec)

	scaleFlo := ring.NewFloat(scale, prec)
	half := ring.NewFloat(0.5, prec)

	for i := range values {

		xFlo.Mul(scaleFlo, values[i])

		if values[i].Cmp(zero) < 0 {
			xFlo.Sub(xFlo, half)
		} else {
			xFlo.Add(xFlo, half)
		}

		xFlo.Int(xInt)

		for j := range moduli {

			Q := ring.NewUint(moduli[j])

			tmp.Mod(xInt, Q)

			if values[i].Cmp(zero) < 0 {
				tmp.Add(tmp, Q)
			}

			coeffs[j][i] = tmp.Uint64()
		}
	}
}

// Divides x by n^2, returns a float
func scaleDown(coeff *big.Int, n float64) (x float64) {

	x, _ = new(big.Float).SetInt(coeff).Float64()
	x /= n

	return
}

func genBigIntChain(Q []uint64) (bigintChain []*big.Int) {

	bigintChain = make([]*big.Int, len(Q))
	bigintChain[0] = ring.NewUint(Q[0])
	for i := 1; i < len(Q); i++ {
		bigintChain[i] = ring.NewUint(Q[i])
		bigintChain[i].Mul(bigintChain[i], bigintChain[i-1])
	}
	return
}

// GenSwitchkeysRescalingParams generates the parameters for rescaling the switching keys
func GenSwitchkeysRescalingParams(Q, P []uint64) (params []uint64) {

	params = make([]uint64, len(Q))

	PBig := ring.NewUint(1)
	for _, pj := range P {
		PBig.Mul(PBig, ring.NewUint(pj))
	}

	tmp := ring.NewUint(0)

	for i := 0; i < len(Q); i++ {

		params[i] = tmp.Mod(PBig, ring.NewUint(Q[i])).Uint64()
		params[i] = ring.ModExp(params[i], int(Q[i]-2), Q[i])
		params[i] = ring.MForm(params[i], Q[i], ring.BRedParams(Q[i]))
	}

	return
}

func sliceBitReverseInPlaceComplex128(slice []complex128, N int) {

	var bit, j int

	for i := 1; i < N; i++ {

		bit = N >> 1

		for j >= bit {
			j -= bit
			bit >>= 1
		}

		j += bit

		if i < j {
			slice[i], slice[j] = slice[j], slice[i]
		}
	}
}

func sliceBitReverseInPlaceRingComplex(slice []*ring.Complex, N int) {

	var bit, j int

	for i := 1; i < N; i++ {

		bit = N >> 1

		for j >= bit {
			j -= bit
			bit >>= 1
		}

		j += bit

		if i < j {
			slice[i], slice[j] = slice[j], slice[i]
		}
	}
}
