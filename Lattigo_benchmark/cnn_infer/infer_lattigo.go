package cnn_infer

import (
	"bufio"
	"fmt"
	"math"
	"os"
	"strconv"
	"time"

	"github.com/ldsec/lattigo/v2/ckks_fv"
)

func Index(slice []int, value int) int {
	for i, v := range slice {
		if v == value {
			return i 
		}
	}
	return -1 
}

func UnitSlice(slice1 []int, slice2 []int) []int {
	for _, rot := range slice2 {
		if Index(slice1, rot) == -1 {
			slice1 = append(slice1, rot)
		}
	}
	return slice1
}

func ResNetCifar10LattigoSparse(layerNum, startImageID, endImageID int) {

	// approximation boundary setting
	B := 40.0;	// approximation boundary
	
	// plain_dnn := false;
	// approx ReLU setting
	alpha := 13		// precision parameter alpha
	comp_no := 3;		// number of compositions
	deg := []int{15,15,27}		// degrees of component polynomials
	scaledVal := 1.7		// scaled_val: the last scaled value
	// double max_factor = 16;		// max_factor = 1 for comparison operation. max_factor > 1 for max or ReLU function
	tree := []*Tree{}		// structure of polynomial evaluation
	evType := oddbaby

	// generate tree
	for i:=0; i<comp_no; i++ {
		tr := Tree{}
		if evType == oddbaby { 
			tr = upgradeOddbaby( deg[i] )  
		} else if evType == baby {
			tr = upgradeBaby( deg[i] )
		} else { panic("evaluation type is not correct") }
		tree = append(tree, &tr)
		// tr.print();
	}

	var outShare *os.File
	var err error

	// Open the file based on layer number
	basePath := "../../cnn_infer/result/"
	filename := basePath + fmt.Sprintf("resnet%d_cifar10_label_%d_%d", layerNum, startImageID, endImageID)
	outShare, err = os.Create(filename)
	if err != nil {
		fmt.Println("Error opening file layer_num is not correct:", err)
		return
	}
	defer outShare.Close()

	endNum := 0		 
	if layerNum == 20 { 
		endNum = 2 // 0 ~ 2
	} else if layerNum == 32 {
		endNum = 4 // 0 ~ 4
	} else if layerNum == 44 {
		endNum = 6	// 0 ~ 6
	} else if layerNum == 56{
		endNum = 8	// 0 ~ 8
	} else if layerNum == 110{
		endNum = 17	// 0 ~ 17
	} else { panic("layer_num is not correct") }

	fmt.Println(endNum)

	// var plaintext *ckks_fv.Plaintext
	btpParams := ckks_fv.DefaultBootstrapParams[0]
	paramsTemp, err := btpParams.Params()
	if err != nil {
		panic(err)
	}

	fmt.Println()
	fmt.Printf("CKKS parameters: logN = %d, logSlots = %d, h = %d, logQP = %d, levels = %d, scale= 2^%f, sigma = %f \n", paramsTemp.LogN(), paramsTemp.LogSlots(), btpParams.H, paramsTemp.LogQP(), paramsTemp.Levels(), math.Log2(paramsTemp.Scale()), paramsTemp.Sigma())

	// Scheme context and keys
	context, _ := genTestParams(paramsTemp, btpParams.H) // generate testContext struc
	fmt.Println()
	fmt.Println("Generating bootstrapping keys...")
	galStepsVector := []int{0}
	for i:=0;i<context.params.LogN()-1;i++{
		galStepsVector = append(galStepsVector, 1<<i)
	}
	rotations := []int {
		0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,
		// ,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55
		56,
		// ,57,58,59,60,61
		62,63,64,66,84,124,128,132,256,512,959,960,990,991,1008,
		1023,1024,1036,1064,1092,1952,1982,1983,2016,2044,2047,2048,2072,2078,2100,3007,3024,3040,3052,3070,3071,3072,3080,3108,4031,
		4032,4062,4063,4095,4096,5023,5024,5054,5055,5087,5118,5119,5120,6047,6078,6079,6111,6112,6142,6143,6144,7071,7102,7103,7135,
		7166,7167,7168,8095,8126,8127,8159,8190,8191,8192,9149,9183,9184,9213,9215,9216,10173,10207,10208,10237,10239,10240,11197,11231,
		11232,11261,11263,11264,12221,12255,12256,12285,12287,12288,13214,13216,13246,13278,13279,13280,13310,13311,13312,14238,14240,
		14270,14302,14303,14304,14334,14335,15262,15264,15294,15326,15327,15328,15358,15359,15360,16286,16288,16318,16350,16351,16352,
		16382,16383,16384,17311,17375,18335,18399,18432,19359,19423,20383,20447,20480,21405,21406,21437,21469,21470,21471,21501,21504,
		22429,22430,22461,22493,22494,22495,22525,22528,23453,23454,23485,23517,23518,23519,23549,24477,24478,24509,24541,24542,24543,
		24573,24576,25501,25565,25568,25600,26525,26589,26592,26624,27549,27613,27616,27648,28573,28637,28640,28672,29600,29632,29664,
		29696,30624,30656,30688,30720,31648,31680,31712,31743,31744,31774,32636,32640,32644,32672,32702,32704,32706,32735,
		32736,32737,32759,32760,32761,32762,32763,32764,32765,32766,32767,
	}
	galStepsVector = UnitSlice(galStepsVector, rotations)

	rotationsBTP := context.kgen.GenRotationIndexesForBootstrapping(15, btpParams)

	galStepsVector = UnitSlice(galStepsVector, rotationsBTP)

	// sort.Ints(rotations) // rotations indices must under order
	fmt.Println("Repeated free rotations: ", galStepsVector)
	rotkeys := context.kgen.GenRotationKeysForRotations(galStepsVector, true, context.sk)
	evk := ckks_fv.EvaluationKey{Rlk: context.rlk, Rtks: rotkeys}

	// context.params.SetLogSlots(15)
	context.evaluator = ckks_fv.NewCKKSEvaluator(context.params, evk)
	var btp *ckks_fv.Bootstrapper
	// var btp2, btp3	*ckks_fv.Bootstrapper
	// context.params.SetLogSlots(14)
	// btpParams.SetLogSlots(14)
	if btp, err = ckks_fv.NewBootstrapperReal(context.params, btpParams, ckks_fv.BootstrappingKey(evk) ); err != nil {
		fmt.Println("Bootstrapper 1 initial failed")
		panic(err)
	}
	// context.params.SetLogSlots(13)
	// btpParams.SetLogSlots(13)
	// if btp2, err = ckks_fv.NewBootstrapper(context.params, btpParams, ckks_fv.BootstrappingKey(evk) ); err != nil {
	// 	fmt.Println("Bootstrapper 2 initial failed")
	// 	panic(err)
	// }
	// context.params.SetLogSlots(12)
	// btpParams.SetLogSlots(12)
	// if btp3, err = ckks_fv.NewBootstrapper(context.params, btpParams, ckks_fv.BootstrappingKey(evk) ); err != nil {
	// 	fmt.Println("Bootstrapper 3 initial failed")
	// 	panic(err)
	// }
	// fmt.Println("Done initial bootstrapper")
	// context.params.SetLogSlots(14)
	// btpParams.SetLogSlots(14)

	// Generate a random plaintext
	valuesWant := make([]complex128, 1<<12 )
	// fmt.Println("\n\nSlots: ", 1<<12 ) 
	for i:=0; i<1<<12;i++ {
		valuesWant[i] = 0.5
	}

	plaintext := context.encoder.EncodeComplexNew(valuesWant, context.params.LogSlots())

	// Encrypt
	ciphertext1 := context.encryptorPk.EncryptNew(plaintext)

	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("Precision of values vs. ciphertext")
	printDebug(ciphertext1, context, 7,7)

	// Bootstrap the ciphertext (homomorphic re-encryption)
	// It takes a ciphertext at level 0 (if not at level 0, then it will reduce it to level 0)
	// and returns a ciphertext at level MaxLevel - k, where k is the depth of the bootstrapping circuit.
	// CAUTION: the scale of the ciphertext MUST be equal (or very close) to params.Scale
	// To equalize the scale, the function evaluator.SetScale(ciphertext, parameters.Scale) can be used at the expense of one level.
	// Start := time.Now()
	ciphertext1 = btp.BootstrappReal(ciphertext1)
	// elapsed := time.Since(Start) // 
	// fmt.Printf("The bootstrap operation took %v\n", elapsed)
	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("Precision of ciphertext vs. Bootstrapp(ciphertext)")
	printDebug(ciphertext1, context, 7,7)
	// // panic("Done ")

	// //****************************

	// context.params.SetLogSlots(15)
	ResidualModuli := ckks_fv.DefaultBootstrapParams[0].ResidualModuli
	fmt.Println("totalLevel", context.params.QiCount(), "ResidualModuli: ",len( ResidualModuli) )
	context.remainLevel = len( ResidualModuli ) - 1 
	logn := context.params.LogSlots()
	fmt.Println("LogSlots: ", logn)

	Start := time.Now()
	for imageID := startImageID; imageID <=endImageID; imageID++ {

		fmt.Println("image id: ", imageID)
		var output *os.File
		var err error

		// Open the file based on layer number
		basePath := "../../cnn_infer/result/"
		filename := basePath + fmt.Sprintf("resnet%d_cifar10_image_%d.txt", layerNum, imageID)
		output, err = os.Create(filename)
		if err != nil {
			fmt.Println("Error opening file layer_num is not correct:", err)
			return
		}
		defer output.Close() 

		// dir := fmt.Sprintf("resnet%d_new", layerNum)

		cipherPool := make( []*ckks_fv.Ciphertext, 14 )

		var cnn *TensorCipher
		var temp *TensorCipher
		// deep learning parameters and import
		var co, st int
		fh, fw := 3, 3

		initP := 8
		n := 1<<logn
		stage := 0
		epsilon := 0.00001
		var image, linearWeight, linearBias []float64
		var convWeight, bnBias, bnRunningMean, bnRunningVar, bnWeight [][]float64

		fmt.Println("Importing parameters")
		importParametersCifar10(&linearWeight, &linearBias, &convWeight, &bnBias, &bnRunningMean, &bnRunningVar, &bnWeight, layerNum, endNum)

		// Pack images compactly
		var in *os.File
		in, err = os.Open("../../cnn_infer/testFile/test_values.txt")
		if err != nil {
			fmt.Println("Error opening file:", err)
			return
		}
		defer in.Close()

		scanner := bufio.NewScanner(in)
		scanner.Split(bufio.ScanWords)
		// Move to the correct image ID
		for i := 0; i < 32*32*3*imageID && scanner.Scan(); i++ { }
		// Set image values to actual values
		image = make([]float64, n)
		for i := 0; i < 32*32*3 && scanner.Scan(); i++ {
			val, _ := strconv.ParseFloat(scanner.Text(), 64)
			image[i] = val
		}

		// Repeat image
		for i := n/initP; i < n; i++ {
			image[i] = image[ i % (n/initP) ]
		}

		// Divide image by B
		for i := 0; i < n; i++ {
			image[i] /= B // For boundary [-1,1]
		}

		// Open the file with image labels
		labelFile, err := os.Open("../../cnn_infer/testFile/test_label.txt")
		if err != nil {
			fmt.Println("Error opening file:", err)
			return
		}
		defer labelFile.Close()

		// Read the label for the current image ID
		scanner = bufio.NewScanner(labelFile)
		scanner.Split(bufio.ScanWords)
		var imageLabel int
		for i := 0; i <= imageID && scanner.Scan(); i++ {
			imageLabel, _ = strconv.Atoi(scanner.Text())
		}
		fmt.Println("image label:", imageLabel)

		// Generate CIFAR-10 image
		// Assuming logn, 1, 32, 32, 3, 3, B, image
		// are the correct parameters for NewTensor.
		context.params.SetScale( float64(1<<51) )
		cnn, _ = NewTensorCipherWithData(logn, 1, 32, 32, 3, 3, initP, image, context, 51 )
		fmt.Println("Done generating Tensor")
		// layer 0
		fmt.Println("layer 0")

		context.evaluator.DropLevel(cnn.cipher, context.totalLevel-context.remainLevel-3 )
		
		cnn, _ = MultiplexedParallelConvolutionCipher(cnn, 16, 1, fh, fw, convWeight[stage], bnRunningVar[stage], bnWeight[stage], epsilon, context, cipherPool, false)
		
		fmt.Println("after convolution")
		printDebug( cnn.cipher, context, 7 ,7 )

		// // scaling factor ~2^51 -> 2^46
		context.params.SetScale( float64(1<<46) )
		modulus := context.ringQ.Modulus
		fmt.Println("modulus: ", math.Log2(float64(modulus[cnn.cipher.Level()])), "cnn.cipher Scale: ", math.Log2(cnn.cipher.Scale()), "cur Scale: ", math.Log2(context.params.Scale()) )
		scaleChange := context.params.Scale() * float64( modulus[cnn.cipher.Level()] ) / (cnn.cipher.Scale())
		plain := ckks_fv.NewPlaintextCKKS(context.params, cnn.cipher.Level(), scaleChange )
		constOne := make([]complex128, context.params.Slots())
		for i := range constOne {
			constOne[i] = complex(1.0, 0.0)
		}
		context.encoder.EncodeComplexNTT( plain, constOne, context.params.LogSlots() )
		context.evaluator.Mul(cnn.cipher, plain, cnn.cipher)
		context.evaluator.RescaleMany(cnn.cipher, 1, cnn.cipher)
		fmt.Println("after Modulus change")
		printDebug( cnn.cipher, context, 7 ,7 )

		cnn, _ = MultiplexedParallelBatchNormCipher(cnn, bnBias[stage], bnRunningMean[stage], bnRunningVar[stage], bnWeight[stage], epsilon, B, context, false)
		fmt.Println("after BN")
		printDebug( cnn.cipher, context, 7 ,7 )

		cnn, _ = ReLU_cipher(cnn, comp_no, deg, alpha, tree, scaledVal, context, btp)
		fmt.Println("after ReLU_cipher")
		printDebug( cnn.cipher, context, 7 ,7 )


		for j:=0; j<3; j++ {// layer 1_x, 2_x, 3_x
		
			if j==0 { 
				co = 16 
			} else if j==1{
				co = 32 
			} else if j==2 { 
				co = 64
			}
			for k:=0; k<=endNum; k++ 	{// 0 ~ 2/4/6/8/17
				stage = 2*((endNum+1)*j+k)+1 
				fmt.Println("layer: ", stage)

				temp = cnn

				if j>=1 && k==0 {
					st = 2
				} else{
					st = 1
				}
				
				cnn, _ = MultiplexedParallelConvolutionCipher(cnn, co, st, fh, fw, convWeight[stage], bnRunningVar[stage], bnWeight[stage], epsilon, context, cipherPool, false)
				fmt.Println("after convolution layer: ",  stage)
				printDebug( cnn.cipher, context, 7 ,7 )

				cnn, _ = MultiplexedParallelBatchNormCipher(cnn, bnBias[stage], bnRunningMean[stage], bnRunningVar[stage], bnWeight[stage], epsilon, B, context, false)
				fmt.Println("after BN layer: ",  stage)
				printDebug( cnn.cipher, context, 7 ,7 )

				// if j==0 { 
				// 	context.params.SetLogSlots(14)
				// 	btpParams.SetLogSlots(14)
				// 	cnn.logn = 14
				// 	cnn.cipher = btp.Bootstrapp(cnn.cipher) 
				// } else if j == 1 {
				// 	context.params.SetLogSlots(13)
				// 	btpParams.SetLogSlots(13)
				// 	cnn.logn = 13
				// 	cnn.cipher = btp2.Bootstrapp(cnn.cipher) 
				// } else {
				// 	context.params.SetLogSlots(12)
				// 	btpParams.SetLogSlots(12)
				// 	cnn.logn = 12
				// 	cnn.cipher = btp3.Bootstrapp(cnn.cipher) 
				// }
				// cnn.cipher = bootstrapReEnc(cnn.cipher, context)
				// fmt.Println("after BTS reEnc: ",  stage)
				// printDebug( cnn.cipher, context, 7 ,7 )
				cnn.cipher.SetScale(context.params.Scale())
				cnn.cipher = btp.BootstrappReal(cnn.cipher) 
				fmt.Println("after BTS layer: ",  stage)
				printDebug( cnn.cipher, context, 7 ,7 )

				cnn, _ = ReLU_cipher(cnn, comp_no, deg, alpha, tree, scaledVal, context, btp)
				fmt.Println("after ReLU_cipher")
				printDebug( cnn.cipher, context, 7 ,7 )

				stage = 2*((endNum+1)*j+k)+2	
				fmt.Println("layer: ", stage)
				st = 1;
				cnn, _ = MultiplexedParallelConvolutionCipher(cnn, co, st, fh, fw, convWeight[stage], bnRunningVar[stage], bnWeight[stage], epsilon, context, cipherPool, false)
				fmt.Println("after convolution layer: ",  stage)
				printDebug( cnn.cipher, context, 7 ,7 )

				cnn, _ = MultiplexedParallelBatchNormCipher(cnn, bnBias[stage], bnRunningMean[stage], bnRunningVar[stage], bnWeight[stage], epsilon, B, context, false)
				fmt.Println("after BN layer: ",  stage)
				printDebug( cnn.cipher, context, 7 ,7 )

				if j>=1 && k==0 {
					temp, _ = multiplexedParallelDownsamplingCipher( temp, context )
					fmt.Println("after Down sample layer: ",  stage)
					printDebug(temp.cipher, context, 7, 7)	
						
				}
				
				cnn, _ = cnnAddCipher(cnn, temp, context)
				fmt.Println("after cnnAddCipher layer: ",  stage)
				printDebug( cnn.cipher, context, 7, 7 )
				
				// if j==0 { 
				// 	context.params.SetLogSlots(14)
				// 	btpParams.SetLogSlots(14)
				// 	cnn.logn = 14
				// 	cnn.cipher = btp.Bootstrapp(cnn.cipher) 
				// } else if j == 1 {
				// 	context.params.SetLogSlots(13)
				// 	btpParams.SetLogSlots(13)
				// 	cnn.logn = 13
				// 	cnn.cipher = btp2.Bootstrapp(cnn.cipher) 
				// } else {
				// 	context.params.SetLogSlots(12)
				// 	btpParams.SetLogSlots(12)
				// 	cnn.logn = 12
				// 	cnn.cipher = btp3.Bootstrapp(cnn.cipher) 
				// }
				// cnn.cipher = bootstrapReEnc(cnn.cipher, context)
				// fmt.Println("after BTS reEnc: ",  stage)
				// printDebug( cnn.cipher, context, 7 ,7 )
				cnn.cipher.SetScale(context.params.Scale())
				cnn.cipher = btp.BootstrappReal(cnn.cipher) 
				fmt.Println("after BTS layer: ",  stage)
				printDebug( cnn.cipher, context, 7 ,7 )
				cnn, _ = ReLU_cipher(cnn, comp_no, deg, alpha, tree, scaledVal, context, btp)
				fmt.Println("after ReLU_cipher")
				printDebug( cnn.cipher, context, 7 ,7 )
			}

		}

		fmt.Println("layer: ", layerNum-1)

		cnn, _ = averagePoolingCipherScale(cnn, B, context)
		fmt.Println("after averagePoolingCipherScale ")
		printDebug( cnn.cipher, context, 7 ,7 )

		cnn, _ = matrixMultiplicationCipher(cnn, linearWeight, linearBias, 10, 64, context)
		// PrintDebugVec(cnn.vec, 7, 7)
		fmt.Println("after matrixMultiplicationCipher ")
		printDebug( cnn.cipher, context, 7 ,7 )

		rtnVec := context.encoder.DecodeComplex(context.decryptor.DecryptNew(cnn.cipher), context.params.LogSlots())
		
		fmt.Println("( ")
		for i:=0; i<10; i++{
			fmt.Println(rtnVec[i])
		}
		fmt.Println(" )")
		

		label := 0
		maxScore := -100.0
		for i:=0; i<10; i++ {
			if maxScore < real(rtnVec[i]) {
				label = i
				maxScore = real(rtnVec[i])
			}
		}
		fmt.Println("image label: ", imageLabel)
		fmt.Println("inferred label: ", label)
		fmt.Println("max score: ", maxScore)

	}
	elapsed := time.Since(Start) 
	fmt.Printf("The infer operation took %v\n", elapsed)

}

func printDebug(ciphertext *ckks_fv.Ciphertext, context *testParams, start, end int) {

	valuesTest := context.encoder.DecodeComplex( context.decryptor.DecryptNew(ciphertext), context.params.LogSlots() )

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), context.params.LogQLvl(ciphertext.Level()))
	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale()))
	fmt.Printf("Level: %d\n", ciphertext.Level())
	fmt.Print("ValuesTest: [")
	for _, v := range valuesTest[:start] {
		fmt.Printf("%f ", real(v) )
		// if i%(33*33*3)==0 {fmt.Println("\n\n")}
	}
	fmt.Print("... ")
	for _, v := range valuesTest[len(valuesTest)-end:] {
		fmt.Printf("%f ", real(v) )
		// if i%(33*33*3)==0 {fmt.Println("\n\n")}
	}
	fmt.Println("]")
	// precStats := ckks_fv.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, params.LogSlots(), 0)

	// fmt.Println(precStats.String())
	fmt.Println()
}
