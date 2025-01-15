package minimax_comp

import (
	"os"
	"strconv"
	"strings"

	"github.com/ldsec/lattigo/v2/ckks_fv"
)

func minimaxComp(compNo int, deg []int, alpha int, trees []*Tree, scaledVal float64, context *testParams, ct0, ct1 *ckks_fv.Ciphertext) (*ckks_fv.Ciphertext, error) {
	decompCoeff := make([][]float64, compNo)
	scaleVal := make([]float64, compNo)

	var tempCipher, halfCipher, xCipher *ckks_fv.Ciphertext

	scaleVal[0] = 1.0
	for i := 1; i < compNo; i++ {
		scaleVal[i] = 2.0
	}
	scaleVal[compNo-1] = scaledVal

	addr := "../../minimax_comp/result/"
	filePath := addr + "d" + strconv.Itoa(alpha) + ".txt"
	data, err := os.ReadFile(filePath)
	if err != nil {
		panic("Read Approximation sgn Coefficient failed!\n")
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
	// (a-b)
	xCipher = context.evaluator.SubNew(ct0, ct1)
	subCipher := xCipher.CopyNew().Ciphertext()
	// (a+b)
	sumCipher := context.evaluator.AddNew(ct0, ct1)

	complexVec := make([]complex128, len(mHalf))
	for i, value := range mHalf { complexVec[i] = complex(value, 0.0) }
	plainHalf := ckks_fv.NewPlaintextCKKS( context.params, sumCipher.Level(), sumCipher.Scale() )
	context.encoder.EncodeComplex( plainHalf, complexVec, context.params.LogSlots())
	halfCipher = context.encryptorPk.EncryptNew( plainHalf )
	// (a+b)/2
	tempCipher = context.evaluator.MulRelinNew(sumCipher, halfCipher )
	context.evaluator.RescaleMany(tempCipher, 1, tempCipher)

	for i := 0; i < compNo; i++ {
		// fmt.Println("in poly evaluate\n",i)
		xCipher = evalPolynomialIntegrateCipher( xCipher, deg[i], decompCoeff[i], trees[i], context)
		// printDebug(xCipher, context, 7, 7)
	}

	// (a+b)+xsgn(x) / 2 from sgn(x)/2
	context.evaluator.MulRelin(subCipher, xCipher, xCipher)
	context.evaluator.RescaleMany(xCipher, 1, xCipher)
	context.evaluator.Add(xCipher, tempCipher, tempCipher)

	return tempCipher, nil
}

func minimaxSort(compNo int, deg []int, alpha int, trees []*Tree, scaledVal float64, context *testParams, ct0, ct1 *ckks_fv.Ciphertext) (*ckks_fv.Ciphertext, *ckks_fv.Ciphertext) {
	decompCoeff := make([][]float64, compNo)
	scaleVal := make([]float64, compNo)

	var tempCipher, halfCipher, xCipher *ckks_fv.Ciphertext

	scaleVal[0] = 1.0
	for i := 1; i < compNo; i++ {
		scaleVal[i] = 2.0
	}
	scaleVal[compNo-1] = scaledVal

	addr := "../../minimax_comp/result/"
	filePath := addr + "d" + strconv.Itoa(alpha) + ".txt"
	data, err := os.ReadFile(filePath)
	if err != nil {
		panic("Read Approximation sgn Coefficient failed!\n")
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
	// (a-b)
	xCipher = context.evaluator.SubNew(ct0, ct1)
	subCipher := xCipher.CopyNew().Ciphertext()
	// (a+b)
	sumCipher := context.evaluator.AddNew(ct0, ct1)

	complexVec := make([]complex128, len(mHalf))
	for i, value := range mHalf { complexVec[i] = complex(value, 0.0) }
	plainHalf := ckks_fv.NewPlaintextCKKS( context.params, sumCipher.Level(), sumCipher.Scale() )
	context.encoder.EncodeComplex( plainHalf, complexVec, context.params.LogSlots())
	halfCipher = context.encryptorPk.EncryptNew( plainHalf )
	// (a+b)/2
	tempCipher = context.evaluator.MulRelinNew(sumCipher, halfCipher )
	context.evaluator.RescaleMany(tempCipher, 1, tempCipher)

	for i := 0; i < compNo; i++ {
		// fmt.Println("in poly evaluate\n",i)
		xCipher = evalPolynomialIntegrateCipher( xCipher, deg[i], decompCoeff[i], trees[i], context)
		// printDebug(xCipher, context, 7, 7)
	}

	// (a+b)+xsgn(x) / 2 from sgn(x)/2
	context.evaluator.MulRelin(subCipher, xCipher, xCipher)
	context.evaluator.RescaleMany(xCipher, 1, xCipher)
	context.evaluator.Add(xCipher, tempCipher, tempCipher)

	return tempCipher, context.evaluator.SubNew(sumCipher, tempCipher)
}


