package cnn_infer

import (
	"bufio"
	"crypto/rand"
	"fmt"
	"math"
	"os"
	"strconv"
	"time"

	"github.com/ldsec/lattigo/v2/ckks_fv"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"golang.org/x/crypto/sha3"
)


func ResNetCifar10RtfRubatoInfer(layerNum, startImageID, endImageID int, name string, rubatoParam int) {

	// approximation boundary setting
	B := 40.0	// approximation boundary
	init_p := 8
	// plain_dnn := false;
	// approx ReLU setting
	alpha := 13		// precision parameter alpha
	comp_no := 3		// number of compositions
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
	btpParams := ckks_fv.DefaultBootstrapParams[0] // generate btp parameters 
	hbtpParams := btpParams.HalfbtpParams() // get hbtp parameter from btp
	paramsTemp, err := btpParams.Params() // get parameters from btp parameters
	if err != nil {
		panic(err)
	}

	// rubato parameters in RtF
	// Rubato parameter
	blocksize := ckks_fv.RubatoParams[rubatoParam].Blocksize
	outputsize := blocksize - 4
	numRound := ckks_fv.RubatoParams[rubatoParam].NumRound
	plainModulus := ckks_fv.RubatoParams[rubatoParam].PlainModulus
	sigma := ckks_fv.RubatoParams[rubatoParam].Sigma


	var fvEncoder ckks_fv.MFVEncoder
	var fvEncryptor ckks_fv.MFVEncryptor
	var fvEvaluator ckks_fv.MFVEvaluator
	paramsTemp.SetPlainModulus(plainModulus)
	paramsTemp.SetLogFVSlots(paramsTemp.LogN()) //using full coeffs
	messageScaling := float64(paramsTemp.PlainModulus()) / hbtpParams.MessageRatio

	fmt.Println()
	fmt.Printf("CKKS parameters: logN = %d, logSlots = %d, h = %d, logQP = %d, levels = %d, scale= 2^%f, sigma = %f, message scaling(log2) = %f \n", paramsTemp.LogN(), paramsTemp.LogSlots(), btpParams.H, paramsTemp.LogQP(), paramsTemp.Levels(), math.Log2(paramsTemp.Scale()), paramsTemp.Sigma(), math.Log2(messageScaling))
	
	rubatoModDown := ckks_fv.RubatoModDownParams[rubatoParam].CipherModDown
	stcModDown := ckks_fv.RubatoModDownParams[rubatoParam].StCModDown
	// Scheme context and keys
	context, _ := genTestParams(paramsTemp, btpParams.H) // generate testContext struc
	fmt.Println()
	fmt.Println("Generating bootstrapping keys...")

	fvEncoder = ckks_fv.NewMFVEncoder(context.params)
	fvEncryptor = ckks_fv.NewMFVEncryptorFromPk(context.params, context.pk)

	// Generating half-bootstrapping keys
	rotationsHalfBoot := context.kgen.GenRotationIndexesForHalfBoot(context.params.LogSlots(), hbtpParams)
	pDcds := fvEncoder.GenSlotToCoeffMatFV(2) // radix 2
	rotationsStC := context.kgen.GenRotationIndexesForSlotsToCoeffsMat(pDcds)
	rotationsRtF := append(rotationsHalfBoot, rotationsStC...)

 	//************************************* 
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
		
		// 4096,8192,12288,16384,20480,24576,28672,
		// 3072,6144,9216,12288,15360,18432,21504,24576,27648,  	//rotate2top 32*32*3*i
	}
	galStepsVector = UnitSlice(galStepsVector, rotations)
	rotationsBTP := context.kgen.GenRotationIndexesForBootstrapping(15, btpParams)
	galStepsVector = UnitSlice(galStepsVector, rotationsBTP)
	galStepsVector = UnitSlice(galStepsVector, rotationsRtF) //get RtF rotation key
	// sort.Ints(rotations) // rotations indices must under order
	fmt.Println("Repeated free rotations: ", galStepsVector)
	rotkeys := context.kgen.GenRotationKeysForRotations(galStepsVector, true, context.sk)

	evk := ckks_fv.EvaluationKey{Rlk: context.rlk, Rtks: rotkeys}

	fvEvaluator = ckks_fv.NewMFVEvaluator(context.params, evk, pDcds)
	context.evaluator = ckks_fv.NewCKKSEvaluator(context.params, evk)
	var btp *ckks_fv.Bootstrapper
	var hbtp *ckks_fv.HalfBootstrapper

	if btp, err = ckks_fv.NewBootstrapperReal(context.params, btpParams, ckks_fv.BootstrappingKey(evk) ); err != nil {
		fmt.Println("Bootstrapper 1 initial failed")
		panic(err)
	}
	if hbtp, err = ckks_fv.NewHalfBootstrapper(context.params, hbtpParams, ckks_fv.BootstrappingKey(evk) ); err != nil {
		panic(err)
	}

	// RtF transciphering framework***********************
	var plainCKKSRingTs []*ckks_fv.PlaintextRingT
	var plaintexts []*ckks_fv.Plaintext
	var rubato ckks_fv.MFVRubato

	var data [][]float64
	var nonces [][]byte
	var counter []byte
	var key []uint64
	var keystream [][]uint64
	var fvKeystreams []*ckks_fv.Ciphertext

	coeffs := make([][]float64, outputsize)
	for s := 0; s < outputsize; s++ {
		coeffs[s] = make([]float64, context.params.N())
	}

	key = make([]uint64, blocksize)
	for i := 0; i < blocksize; i++ {
		key[i] = uint64(i + 1) // Use (1, ..., 16) for testing
	}
	data = make([][]float64, outputsize)

	// Load images parameters from files*********************************
	var imagesCipher []*ckks_fv.Ciphertext
	formattedImage := make([][][]float64, outputsize) // row 8, column 20 store total 80 images
	
	for i := range formattedImage {
		formattedImage[i] = make([][]float64, 20) //20 images each row 
		for j := range formattedImage[i] {
			formattedImage[i][j] = make([]float64, 32*32*3)
		}
	}
	for s := 0; s < outputsize; s++ {
		data[s] = make([]float64, context.params.N())
	}

	if true {
		image_num := endImageID-startImageID+1 //total images
		
		if image_num > outputsize * 20 {panic("images loaded is oversize, that is outputsize * 20 for fullCoeffs")}

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

			image := make([]float64, 32*32*3)

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
			for i := 0; i < 32*32*3 && scanner.Scan(); i++ {
				val, _ := strconv.ParseFloat(scanner.Text(), 64)
				image[i] = val
			}
			// Divide image by B
			for i := 0; i < len(image); i++ {
				image[i] /= B 	// For boundary [-1,1]
			}
			//debug
			if  imageID - startImageID == 0 {
				fmt.Println("plaintext data image0 = ")
				for i:=0;i<20;i++ {
					fmt.Println("(", image[i], "), ")
					if i%20 == 0 {fmt.Println("")}
				}
			}
			var cur_col, cur_row int
			cur_col = (imageID-startImageID) % 20
			cur_row = (imageID-startImageID) / 20

			copy(formattedImage[cur_row][cur_col], image)	
		}
		// first part data embadding 
		for i:=0; i<10*(32*32*3); i++{
			cur_image_col := i / (32*32*3)
			for row:=0; row<outputsize; row++ {
				data[row][i] = formattedImage[row][cur_image_col][i%(32*32*3)]
			}
		}
		// second part of data embadding
		for i:=context.params.Slots(); i<context.params.Slots() + 10*(32*32*3); i++{
			cur_image_col := (i-context.params.Slots()) / (32*32*3) + 10
			for row:=0; row<outputsize; row++ {
				data[row][i] = formattedImage[row][cur_image_col][i%(32*32*3)]
			}
		}	

	} else {
		fmt.Println("Do not enable transciphering!\n ")
	}
	// Load images parameters from files*********************************
	fmt.Println("fullSlots, function")
	
	nonces = make([][]byte, context.params.N())
	for i := 0; i < context.params.N(); i++ {
		nonces[i] = make([]byte, 8)
		rand.Read(nonces[i])
	}
	counter = make([]byte, 8)
	rand.Read(counter)

	// Get keystream, focus on the coefficient format of the encoding 
	keystream = make([][]uint64, context.params.N())
	for i := 0; i < context.params.N(); i++ {
		keystream[i] = plainRubato(blocksize, numRound, nonces[i], counter, key, plainModulus, sigma)
	}

	for s := 0; s < outputsize; s++ {
		for i := 0; i < context.params.N()/2; i++ {
			j := utils.BitReverse64(uint64(i), uint64(context.params.LogN()-1))
			coeffs[s][j] = data[s][i]
			coeffs[s][j+uint64(context.params.N()/2)] = data[s][i+context.params.N()/2]
		}
	}

	// Encode plaintext	
	plainCKKSRingTs = make([]*ckks_fv.PlaintextRingT, 16)
	for s := 0; s < outputsize; s++ {
		plainCKKSRingTs[s] = context.encoder.EncodeCoeffsRingTNew(coeffs[s], messageScaling) //encode plain into coeff
		poly := plainCKKSRingTs[s].Value()[0] // RNS of RingTs has only one moduli
		for i := 0; i < context.params.N(); i++ {
			j := utils.BitReverse64(uint64(i), uint64(context.params.LogN()))
			poly.Coeffs[0][j] = (poly.Coeffs[0][j] + keystream[i][s]) % context.params.PlainModulus()
		}
	}

	plaintexts = make([]*ckks_fv.Plaintext, outputsize) // plaintext used to store plainCKKSRingTs obtained from symmetric

	for s := 0; s < outputsize; s++ {
		plaintexts[s] = ckks_fv.NewPlaintextFVLvl(context.params, 0) // generate a level one FV plaintext
		fvEncoder.FVScaleUp(plainCKKSRingTs[s], plaintexts[s]) // transforms a PlaintextRingT (R_t) into a Plaintext (R_q) by scaling up the coefficient by Q/t.
	}

	rubato = ckks_fv.NewMFVRubato(rubatoParam, context.params, fvEncoder, fvEncryptor, fvEvaluator, rubatoModDown[0])
	kCt := rubato.EncKey(key)

	tc := time.Now()
	// FV Keystream
	fvKeystreams = rubato.Crypt(nonces, counter, kCt, rubatoModDown) // homomorphically generate keystream encrypted under FV
	for i := 0; i < 1; i++ {
		fvKeystreams[i] = fvEvaluator.SlotsToCoeffs(fvKeystreams[i], stcModDown) // FV S2C to get Coeff format keystream
		fvEvaluator.ModSwitchMany(fvKeystreams[i], fvKeystreams[i], fvKeystreams[i].Level())
	}

	for i := 1; i < outputsize; i++ {
		fvKeystreams[i] = fvEvaluator.SlotsToCoeffs(fvKeystreams[i], stcModDown)
		fvEvaluator.ModSwitchMany(fvKeystreams[i], fvKeystreams[i], fvKeystreams[i].Level())
	}

	var ctBoot, ctBoot2 *ckks_fv.Ciphertext

	// Encrypt and mod switch to the lowest level
	ciphertext := ckks_fv.NewCiphertextFVLvl(context.params, 1, 0)
	ciphertext.Value()[0] = plaintexts[0].Value()[0].CopyNew()
	fvEvaluator.Sub(ciphertext, fvKeystreams[0], ciphertext)
	fvEvaluator.TransformToNTT(ciphertext, ciphertext)

	ciphertext.SetScale(math.Exp2(math.Round(math.Log2(float64(context.params.Qi()[0]) / float64(context.params.PlainModulus()) * messageScaling))))

	// Half-Bootstrap the ciphertext (homomorphic evaluation of ModRaise -> SubSum -> CtS -> EvalMod)
	// It takes a ciphertext at level 0 (if not at level 0, then it will reduce it to level 0)
	// and returns a ciphertext at level MaxLevel - k, where k is the depth of the bootstrapping circuit.
	// Difference from the bootstrapping is that the last StC is missing.
	// CAUTION: the scale of the ciphertext MUST be equal (or very close) to params.Scale
	// To equalize the scale, the function evaluator.SetScale(ciphertext, parameters.Scale) can be used at the expense of one level.
	ctBoot, ctBoot2 = hbtp.HalfBoot(ciphertext, false)
	endTc := time.Since(tc)
	fmt.Println("Time transcipher : %v", endTc)
	fmt.Println("after HalfBoot: ")
	printDebug(ctBoot, context, 7, 7)
	printDebug(ctBoot2, context, 7, 7)
	
	for i:=0; i<1; i++{
		for pos:=0; pos<10; pos++{
			rtnCipher := ProcessTranscipherResult(ctBoot, init_p, 32*32*3, pos, context)
			imagesCipher = append(imagesCipher, rtnCipher)
		}
	}
	printDebug(imagesCipher[0], context, 7, 7)


	// valuesWant := make([]complex128, context.params.Slots())
	// for i := 0; i < context.params.Slots(); i++ {
	// 	valuesWant[i] = complex(data[0][i], 0)
	// }
	// context.evaluator.Mul(ctBoot, ctBoot, ctBoot)
	// context.evaluator.RescaleMany(ctBoot, 1, ctBoot)
	// fmt.Println("Precision of HalfBoot(ciphertext)")
	// printDebug(ctBoot, context, 7, 7)

	// RtF transciphering framework***********************
	// panic("Done RtF framework")
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
		// context.params.SetScale( float64(1<<51) )
		// cnn, _ = NewTensorCipherWithData(logn, 1, 32, 32, 3, 3, initP, image, context, 51 )
		cnn, _ = NewTensorCipherWithCT(logn, 1, 32, 32, 3, 3, initP, imagesCipher[imageID] )

		fmt.Println("Done generating Tensor")
		// layer 0
		fmt.Println("layer 0")

		for cnn.cipher.Level() > 2 {
			context.evaluator.DropLevel(cnn.cipher, 1)
		}
		
		cnn, _ = MultiplexedParallelConvolutionCipher(cnn, 16, 1, fh, fw, convWeight[stage], bnRunningVar[stage], bnWeight[stage], epsilon, context, cipherPool, false)
		
		fmt.Println("after convolution")
		printDebug( cnn.cipher, context, 7 ,7 )

		cnn, _ = MultiplexedParallelBatchNormCipher(cnn, bnBias[stage], bnRunningMean[stage], bnRunningVar[stage], bnWeight[stage], epsilon, B, context, false)
		fmt.Println("after BN")
		printDebug( cnn.cipher, context, 7 ,7 )

		cnn.cipher.SetScale(context.params.Scale())
		cnn.cipher = btp.BootstrappReal(cnn.cipher)
		fmt.Println("after btp real, stage: ", stage)
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


func plainRubato(blocksize int, numRound int, nonce []byte, counter []byte, key []uint64, plainModulus uint64, sigma float64) (state []uint64) {
	xof := sha3.NewShake256()
	xof.Write(nonce)
	xof.Write(counter)
	state = make([]uint64, blocksize)

	prng, err := utils.NewPRNG()
	if err != nil {
		panic(err)
	}
	gaussianSampler := ring.NewGaussianSampler(prng)

	rks := make([][]uint64, numRound+1)

	for r := 0; r <= numRound; r++ {
		rks[r] = make([]uint64, blocksize)
		for i := 0; i < blocksize; i++ {
			rks[r][i] = ckks_fv.SampleZqx(xof, plainModulus) * key[i] % plainModulus
		}
	}

	for i := 0; i < blocksize; i++ {
		state[i] = uint64(i + 1)
	}

	// Initial AddRoundKey
	for i := 0; i < blocksize; i++ {
		state[i] = (state[i] + rks[0][i]) % plainModulus
	}

	// Round Functions
	for r := 1; r < numRound; r++ {
		rubatoLinearLayer(state, plainModulus)
		rubatoFeistel(state, plainModulus)
		for i := 0; i < blocksize; i++ {
			state[i] = (state[i] + rks[r][i]) % plainModulus
		}
	}

	// Finalization
	rubatoLinearLayer(state, plainModulus)
	rubatoFeistel(state, plainModulus)
	rubatoLinearLayer(state, plainModulus)
	if sigma > 0 {
		rubatoAddGaussianNoise(state, plainModulus, gaussianSampler, sigma)
	}
	for i := 0; i < blocksize; i++ {
		state[i] = (state[i] + rks[numRound][i]) % plainModulus
	}
	state = state[0 : blocksize-4]

	return
}

func rubatoLinearLayer(state []uint64, plainModulus uint64) {
	blocksize := len(state)
	buf := make([]uint64, blocksize)

	if blocksize == 16 {
		// MixColumns
		for row := 0; row < 4; row++ {
			for col := 0; col < 4; col++ {
				buf[row*4+col] = 2 * state[row*4+col]
				buf[row*4+col] += 3 * state[((row+1)%4)*4+col]
				buf[row*4+col] += state[((row+2)%4)*4+col]
				buf[row*4+col] += state[((row+3)%4)*4+col]
				buf[row*4+col] %= plainModulus
			}
		}
		// MixRows
		for row := 0; row < 4; row++ {
			for col := 0; col < 4; col++ {
				state[row*4+col] = 2 * buf[row*4+col]
				state[row*4+col] += 3 * buf[row*4+(col+1)%4]
				state[row*4+col] += buf[row*4+(col+2)%4]
				state[row*4+col] += buf[row*4+(col+3)%4]
				state[row*4+col] %= plainModulus
			}
		}
	} else if blocksize == 36 {
		// MixColumns
		for row := 0; row < 6; row++ {
			for col := 0; col < 6; col++ {
				buf[row*6+col] = 4 * state[row*6+col]
				buf[row*6+col] += 2 * state[((row+1)%6)*6+col]
				buf[row*6+col] += 4 * state[((row+2)%6)*6+col]
				buf[row*6+col] += 3 * state[((row+3)%6)*6+col]
				buf[row*6+col] += state[((row+4)%6)*6+col]
				buf[row*6+col] += state[((row+5)%6)*6+col]
				buf[row*6+col] %= plainModulus
			}
		}
		// MixRows
		for row := 0; row < 6; row++ {
			for col := 0; col < 6; col++ {
				state[row*6+col] = 4 * buf[row*6+col]
				state[row*6+col] += 2 * buf[row*6+(col+1)%6]
				state[row*6+col] += 4 * buf[row*6+(col+2)%6]
				state[row*6+col] += 3 * buf[row*6+(col+3)%6]
				state[row*6+col] += buf[row*6+(col+4)%6]
				state[row*6+col] += buf[row*6+(col+5)%6]
				state[row*6+col] %= plainModulus
			}
		}
	} else if blocksize == 64 {
		// MixColumns
		for row := 0; row < 8; row++ {
			for col := 0; col < 8; col++ {
				buf[row*8+col] = 5 * state[row*8+col]
				buf[row*8+col] += 3 * state[((row+1)%8)*8+col]
				buf[row*8+col] += 4 * state[((row+2)%8)*8+col]
				buf[row*8+col] += 3 * state[((row+3)%8)*8+col]
				buf[row*8+col] += 6 * state[((row+4)%8)*8+col]
				buf[row*8+col] += 2 * state[((row+5)%8)*8+col]
				buf[row*8+col] += state[((row+6)%8)*8+col]
				buf[row*8+col] += state[((row+7)%8)*8+col]
				buf[row*8+col] %= plainModulus
			}
		}
		// MixRows
		for row := 0; row < 8; row++ {
			for col := 0; col < 8; col++ {
				state[row*8+col] = 5 * buf[row*8+col]
				state[row*8+col] += 3 * buf[row*8+(col+1)%8]
				state[row*8+col] += 4 * buf[row*8+(col+2)%8]
				state[row*8+col] += 3 * buf[row*8+(col+3)%8]
				state[row*8+col] += 6 * buf[row*8+(col+4)%8]
				state[row*8+col] += 2 * buf[row*8+(col+5)%8]
				state[row*8+col] += buf[row*8+(col+6)%8]
				state[row*8+col] += buf[row*8+(col+7)%8]
				state[row*8+col] %= plainModulus
			}
		}
	} else {
		panic("Invalid blocksize")
	}
}

func rubatoFeistel(state []uint64, plainModulus uint64) {
	blocksize := len(state)
	buf := make([]uint64, blocksize)

	for i := 0; i < blocksize; i++ {
		buf[i] = state[i]
	}

	for i := 1; i < blocksize; i++ {
		state[i] = (buf[i] + buf[i-1]*buf[i-1]) % plainModulus
	}
}

func rubatoAddGaussianNoise(state []uint64, plainModulus uint64, gaussianSampler *ring.GaussianSampler, sigma float64) {
	bound := int(6 * sigma)
	gaussianSampler.AGN(state, plainModulus, sigma, bound)
}
