package cnn_infer

import (
	"os"
	"strconv"
	"strings"

	"github.com/ldsec/lattigo/v2/ckks_fv"
)

func minimax_ReLU_cipher(compNo int, deg []int, alpha int, trees []*Tree, scaledVal float64, context *testParams, cipherIn *ckks_fv.Ciphertext, btp *ckks_fv.Bootstrapper) (*ckks_fv.Ciphertext, error) {
	decompCoeff := make([][]float64, compNo)
	scaleVal := make([]float64, compNo)

	var tempCipher, halfCipher, xCipher *ckks_fv.Ciphertext

	scaleVal[0] = 1.0
	for i := 1; i < compNo; i++ {
		scaleVal[i] = 2.0
	}
	scaleVal[compNo-1] = scaledVal

	addr := "../../cnn_infer/result/"
	filePath := addr + "d" + strconv.Itoa(alpha) + ".txt"
	data, err := os.ReadFile(filePath)
	if err != nil {
		panic("Read Approximation ReLU Coefficient failed!\n")
	}

	coeffStrings := strings.Fields(string(data))
	coeffIdx := 0
	for i := range decompCoeff {
		decompCoeff[i] = make([]float64, 0)
		for j := 0; j < coeffNumber(deg[i], trees[i]); j++ {
			coeff, convErr := strconv.ParseFloat(coeffStrings[coeffIdx], 64)
			if convErr != nil {
				panic("Read coefficient from file failed!\n")
			}
			decompCoeff[i] = append(decompCoeff[i], coeff)
			coeffIdx++
		}
	}

	// scale coefficients properly so that unnecessary level consumptions do not occur
	for i := 0; i < compNo-1; i++ {
		for j := range decompCoeff[i] {
			decompCoeff[i][j] /= scaleVal[i+1]
		}
	}
	for j := range decompCoeff[compNo-1] {
		decompCoeff[compNo-1][j] *= 0.5
	}

	n := context.params.Slots()

	mHalf := make([]float64, n)
	for i := range mHalf {
		mHalf[i] = 0.5
	}
	// generation of half ciphertext
	// long n = cipher_in.poly_modulus_degree()/2;
	xCipher = cipherIn

	for i := 0; i < compNo-1; i++ {
		// fmt.Println("in poly evaluate\n",i)
		xCipher = evalPolynomialIntegrateCipher( xCipher, deg[i], decompCoeff[i], trees[i], context)
		// printDebug(xCipher, context, 7, 7)
	}
	xCipher.SetScale(context.params.Scale())
	xCipher = btp.BootstrappReal(xCipher)
	printDebug(xCipher, context, 7, 7)

	xCipher = evalPolynomialIntegrateCipher( xCipher, deg[compNo-1], decompCoeff[compNo-1], trees[compNo-1], context)

	// x(1+sgn(x))/2 from sgn(x)/2
	complexVec := make([]complex128, len(mHalf))
	for i, value := range mHalf { complexVec[i] = complex(value, 0.0) }

	plainHalf := ckks_fv.NewPlaintextCKKS( context.params, xCipher.Level(), xCipher.Scale() )
	context.encoder.EncodeComplex( plainHalf, complexVec, context.params.LogSlots())
	halfCipher = context.encryptorPk.EncryptNew( plainHalf )
	tempCipher = context.evaluator.AddNew(xCipher, halfCipher )
	context.evaluator.MulRelin(tempCipher, cipherIn, tempCipher)
	context.evaluator.RescaleMany(tempCipher, 1, tempCipher)
	return tempCipher, nil
}