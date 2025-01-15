package main

import (
	"github.com/ldsec/lattigo/v2/minimax_comp"
)



func main() {
	//     //2 Example run: import CIFAR-10 parameters for ResNet-20
	// 	// cnn_infer.ResNetCifar10Sparse(20, 0, 0)
		// cnn_infer.ResNetCifar10LattigoSparse(20, 0, 0)
		minimax_comp.MinimaxCompHera("80f", 5, 0, 2, true)
		// minimax_comp.MinimaxCompRubato( "Param0", ckks_fv.RUBATO128S)
		// func ResNetCifar10RtfHeraInfer(layerNum, startImageID, endImageID int, name string, numRound int, paramIndex int, radix int, fullCoeffs bool) {

}
	