package main

import (
	"fmt"
	"math"

	"github.com/ldsec/lattigo/v2/ckks_fv"
	"github.com/ldsec/lattigo/v2/cnn_infer"
)

func main() {
//     //2 Example run: import CIFAR-10 parameters for ResNet-20
// 	// cnn_infer.ResNetCifar10Sparse(20, 0, 0)
	// cnn_infer.ResNetCifar10LattigoSparse(20, 0, 0)
	// cnn_infer.ResNetCifar10RtfHeraInfer(20, 0, 0, "80f", 4, 0, 2, true)
	cnn_infer.ResNetCifar10RtfRubatoInfer(20, 0, 0, "Param0", ckks_fv.RUBATO128S)
	// func ResNetCifar10RtfHeraInfer(layerNum, startImageID, endImageID int, name string, numRound int, paramIndex int, radix int, fullCoeffs bool) {


}

func printDebug(params *ckks_fv.Parameters, ciphertext *ckks_fv.Ciphertext, valuesWant []complex128, decryptor ckks_fv.CKKSDecryptor, encoder ckks_fv.CKKSEncoder) (valuesTest []complex128) {

	valuesTest = encoder.DecodeComplex(decryptor.DecryptNew(ciphertext), params.LogSlots())

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale()))
	fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])

	precStats := ckks_fv.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, params.LogSlots(), 0)

	fmt.Println(precStats.String())
	fmt.Println()

	return
}
