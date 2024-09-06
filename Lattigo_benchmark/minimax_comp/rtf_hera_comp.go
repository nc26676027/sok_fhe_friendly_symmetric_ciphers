package minimax_comp

import (
	"crypto/rand"
	"fmt"
	"math"
	"time"

	"github.com/ldsec/lattigo/v2/ckks_fv"
	"github.com/ldsec/lattigo/v2/utils"
	"golang.org/x/crypto/sha3"
)


func MinimaxCompHera( name string, numRound int, paramIndex int, radix int, fullCoeffs bool) {
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

	var err error

	// var plaintext *ckks_fv.Plaintext
	btpParams := ckks_fv.DefaultBootstrapParams[0] // generate btp parameters 
	hbtpParams := btpParams.HalfbtpParams() // get hbtp parameter from btp
	paramsTemp, err := btpParams.Params() // get parameters from btp parameters
	if err != nil {
		panic(err)
	}
	messageScaling := float64(paramsTemp.PlainModulus()) / hbtpParams.MessageRatio

	fmt.Println()
	fmt.Printf("CKKS parameters: logN = %d, logSlots = %d, h = %d, logQP = %d, levels = %d, scale= 2^%f, sigma = %f, message scaling(log2) = %f \n", paramsTemp.LogN(), paramsTemp.LogSlots(), btpParams.H, paramsTemp.LogQP(), paramsTemp.Levels(), math.Log2(paramsTemp.Scale()), paramsTemp.Sigma(), math.Log2(messageScaling))
	
	// HERA parameters in RtF
	var heraModDown, stcModDown []int
	if numRound == 4 {
		heraModDown = ckks_fv.HeraModDownParams80[paramIndex].CipherModDown
		stcModDown = ckks_fv.HeraModDownParams80[paramIndex].StCModDown
	} else {
		heraModDown = ckks_fv.HeraModDownParams128[paramIndex].CipherModDown
		stcModDown = ckks_fv.HeraModDownParams128[paramIndex].StCModDown
	}

	// Scheme context and keys
	context, _ := genTestParams(paramsTemp, btpParams.H) // generate testContext struc
	fmt.Println()
	fmt.Println("Generating bootstrapping keys...")
	var fvEncoder ckks_fv.MFVEncoder
	var fvEncryptor ckks_fv.MFVEncryptor
	var fvEvaluator ckks_fv.MFVEvaluator
	// fullCoeffs denotes whether full coefficients are used for data encoding
	if fullCoeffs {
		context.params.SetLogFVSlots(context.params.LogN())
	} else {
		context.params.SetLogFVSlots(context.params.LogSlots())
	}
	fvEncoder = ckks_fv.NewMFVEncoder(context.params)
	fvEncryptor = ckks_fv.NewMFVEncryptorFromPk(context.params, context.pk)

	// Generating half-bootstrapping keys
	rotationsHalfBoot := context.kgen.GenRotationIndexesForHalfBoot(context.params.LogSlots(), hbtpParams)
	pDcds := fvEncoder.GenSlotToCoeffMatFV(radix)
	rotationsStC := context.kgen.GenRotationIndexesForSlotsToCoeffsMat(pDcds)

	rotationsRtF := append(rotationsHalfBoot, rotationsStC...)
	if !fullCoeffs {
		rotationsRtF = append(rotationsRtF, context.params.Slots()/2)
	}

 	//************************************* 
	galStepsVector := []int{0}
	for i:=0;i<context.params.LogN()-1;i++{
		galStepsVector = append(galStepsVector, 1<<i)
	}
	rotationsBTP := context.kgen.GenRotationIndexesForBootstrapping(15, btpParams)
	galStepsVector = UnitSlice(galStepsVector, rotationsBTP)

	galStepsVector = UnitSlice(galStepsVector, rotationsRtF) //get RtF rotation key
	// sort.Ints(rotations) // rotations indices must under order
	fmt.Println("Repeated free rotations: ", galStepsVector)
	rotkeys := context.kgen.GenRotationKeysForRotations(galStepsVector, true, context.sk)

	evk := ckks_fv.EvaluationKey{Rlk: context.rlk, Rtks: rotkeys}

	fvEvaluator = ckks_fv.NewMFVEvaluator(context.params, evk, pDcds)
	context.evaluator = ckks_fv.NewCKKSEvaluator(context.params, evk)
	var hbtp *ckks_fv.HalfBootstrapper

	// if btp, err = ckks_fv.NewBootstrapperReal(context.params, btpParams, ckks_fv.BootstrappingKey(evk) ); err != nil {
	// 	fmt.Println("Bootstrapper 1 initial failed")
	// 	panic(err)
	// }

	if hbtp, err = ckks_fv.NewHalfBootstrapper(context.params, hbtpParams, ckks_fv.BootstrappingKey(evk) ); err != nil {
		panic(err)
	}

	// RtF transciphering framework***********************
	var plainCKKSRingTs []*ckks_fv.PlaintextRingT
	var plaintexts []*ckks_fv.Plaintext
	var hera ckks_fv.MFVHera

	var data [][]float64
	var nonces [][]byte
	var key []uint64
	var keystream [][]uint64
	var fvKeystreams []*ckks_fv.Ciphertext
	data = make([][]float64, 16)
	coeffs := make([][]float64, 16)

	for s := 0; s < 16; s++ { // 16 polynomials to store 16block of Hera, using row packing 
		coeffs[s] = make([]float64, context.params.N())
	}

	key = make([]uint64, 16)
	for i := 0; i < 16; i++ {
		key[i] = uint64(i + 1) // Use (1, ..., 16) for testing
	}
	for s := 0; s < 16; s++ {
		data[s] = make([]float64, context.params.N())
		for i := 0; i < context.params.N()/2; i++ {
			data[s][i] = 1
			data[s][i+context.params.N()/2] = 2
		}
	}
	// Load images parameters from files*********************************

	if fullCoeffs {
		fmt.Println("fullSlots, function")
		
		nonces = make([][]byte, context.params.N())
		for i := 0; i < context.params.N(); i++ {
			nonces[i] = make([]byte, 64)
			rand.Read(nonces[i])
		}

		keystream = make([][]uint64, context.params.N())
		for i := 0; i < context.params.N(); i++ {
			keystream[i] = plainHera(numRound, nonces[i], key, context.params.PlainModulus())
		}

		fmt.Println("total parmas.N block of symmetric enc, keystream block 0: ", keystream[0])

		for s := 0; s < 16; s++ {
			for i := 0; i < context.params.N()/2; i++ {
				j := utils.BitReverse64(uint64(i), uint64(context.params.LogN()-1))
				coeffs[s][j] = data[s][i] // ?? first half data encoded to front part coeff
				coeffs[s][j+uint64(context.params.N()/2)] = data[s][i+context.params.N()/2]// ?? second half data encoded to behand part coeff 
			}
		}

		plainCKKSRingTs = make([]*ckks_fv.PlaintextRingT, 16)
		for s := 0; s < 16; s++ {
			plainCKKSRingTs[s] = context.encoder.EncodeCoeffsRingTNew(coeffs[s], messageScaling) //encode plain into coeff
			poly := plainCKKSRingTs[s].Value()[0] // RNS of RingTs has only one moduli
			for i := 0; i < context.params.N(); i++ {
				j := utils.BitReverse64(uint64(i), uint64(context.params.LogN()))
				poly.Coeffs[0][j] = (poly.Coeffs[0][j] + keystream[i][s]) % context.params.PlainModulus()
			}
		}
	} else {
		fmt.Println("no fullSlots: ")

		nonces = make([][]byte, context.params.Slots())
		for i := 0; i < context.params.Slots(); i++ {
			nonces[i] = make([]byte, 64)
			rand.Read(nonces[i])
		}

		keystream = make([][]uint64, context.params.Slots())
		for i := 0; i < context.params.Slots(); i++ {
			keystream[i] = plainHera(numRound, nonces[i], key, context.params.PlainModulus())
		}

		for s := 0; s < 16; s++ {
			for i := 0; i < context.params.Slots()/2; i++ {
				j := utils.BitReverse64(uint64(i), uint64(context.params.LogN()-1))
				coeffs[s][j] = data[s][i]
				coeffs[s][j+uint64(context.params.N()/2)] = data[s][i+context.params.Slots()/2]
			}
		}

		plainCKKSRingTs = make([]*ckks_fv.PlaintextRingT, 16)
		for s := 0; s < 16; s++ {
			plainCKKSRingTs[s] = context.encoder.EncodeCoeffsRingTNew(coeffs[s], messageScaling)
			poly := plainCKKSRingTs[s].Value()[0]
			for i := 0; i < context.params.Slots(); i++ {
				j := utils.BitReverse64(uint64(i), uint64(context.params.LogN()))
				poly.Coeffs[0][j] = (poly.Coeffs[0][j] + keystream[i][s]) % context.params.PlainModulus()
			}
		}

	}

	plaintexts = make([]*ckks_fv.Plaintext, 16) // plaintext used to store plainCKKSRingTs obtained from symmetric

	for s := 0; s < 16; s++ {
		plaintexts[s] = ckks_fv.NewPlaintextFVLvl(context.params, 0) // generate a level one FV plaintext
		fvEncoder.FVScaleUp(plainCKKSRingTs[s], plaintexts[s]) // transforms a PlaintextRingT (R_t) into a Plaintext (R_q) by scaling up the coefficient by Q/t.
	}

	hera = ckks_fv.NewMFVHera(numRound, context.params, fvEncoder, fvEncryptor, fvEvaluator, heraModDown[0]) // initial hera decryption object
	kCt := hera.EncKey(key) // obtion the FV encrypted symmetric key of hera
	
	// FV Keystream
	fvKeystreams = hera.Crypt(nonces, kCt, heraModDown) // homomorphically generate keystream encrypted under FV
	for i := 0; i < 1; i++ {
		fvKeystreams[i] = fvEvaluator.SlotsToCoeffs(fvKeystreams[i], stcModDown) // FV S2C to get Coeff format keystream
		fvEvaluator.ModSwitchMany(fvKeystreams[i], fvKeystreams[i], fvKeystreams[i].Level())
	}

	for i := 1; i < 16; i++ {
		fvKeystreams[i] = fvEvaluator.SlotsToCoeffs(fvKeystreams[i], stcModDown) // homomorphically generate keystream encrypted under FV
		fvEvaluator.ModSwitchMany(fvKeystreams[i], fvKeystreams[i], fvKeystreams[i].Level()) // FV S2C to get Coeff format keystream
	}

	ctBoot :=make([]*ckks_fv.Ciphertext, 4)

	// Encrypt and mod switch to the lowest level
	ciphertext := ckks_fv.NewCiphertextFVLvl(context.params, 1, 0)
	ciphertext.Value()[0] = plaintexts[0].Value()[0].CopyNew()
	fvEvaluator.Sub(ciphertext, fvKeystreams[0], ciphertext)
	fvEvaluator.TransformToNTT(ciphertext, ciphertext)

	ciphertext.SetScale(math.Exp2(math.Round(math.Log2(float64(context.params.Qi()[0]) / float64(context.params.PlainModulus()) * messageScaling))))	

	ciphertext2 := ckks_fv.NewCiphertextFVLvl(context.params, 1, 0)
	ciphertext2.Value()[1] = plaintexts[1].Value()[0].CopyNew()
	fvEvaluator.Sub(ciphertext2, fvKeystreams[1], ciphertext2)
	fvEvaluator.TransformToNTT(ciphertext2, ciphertext2)

	ciphertext2.SetScale(math.Exp2(math.Round(math.Log2(float64(context.params.Qi()[0]) / float64(context.params.PlainModulus()) * messageScaling))))	

	// Half-Bootstrap the ciphertext (homomorphic evaluation of ModRaise -> SubSum -> CtS -> EvalMod)
	// It takes a ciphertext at level 0 (if not at level 0, then it will reduce it to level 0)
	// and returns a ciphertext at level MaxLevel - k, where k is the depth of the bootstrapping circuit.
	// Difference from the bootstrapping is that the last StC is missing.
	// CAUTION: the scale of the ciphertext MUST be equal (or very close) to params.Scale
	// To equalize the scale, the function evaluator.SetScale(ciphertext, parameters.Scale) can be used at the expense of one level.
	if fullCoeffs {
		ctBoot[0], ctBoot[1] = hbtp.HalfBoot(ciphertext, false)
		// ctBoot[2], ctBoot[3] = hbtp.HalfBoot(ciphertext2, false)
	} else {
		// ctBoot, _ = hbtp.HalfBoot(ciphertext, true)
		panic("half coeff not coded")
	}

	// RtF transciphering framework***********************

	// panic("Done RtF framework")
	ResidualModuli := ckks_fv.DefaultBootstrapParams[0].ResidualModuli
	fmt.Println("totalLevel", context.params.QiCount(), "ResidualModuli: ",len( ResidualModuli) )
	context.remainLevel = len( ResidualModuli ) - 1 
	logn := context.params.LogSlots()
	fmt.Println("LogSlots: ", logn)

	Start := time.Now()
	//stage 1
	ctxt1, _:= minimaxComp(comp_no, deg, alpha, tree, scaledVal, context, ctBoot[0], ctBoot[1])
	elapsed := time.Since(Start) 

	rtnVec := context.encoder.DecodeComplex(context.decryptor.DecryptNew(ctxt1), context.params.LogSlots())
	
	ckks_fv.PrintDebugVecC(rtnVec, 7, 3)

	fmt.Printf("The infer operation took %v\n", elapsed)
}


func plainHera(roundNum int, nonce []byte, key []uint64, plainModulus uint64) (state []uint64) {
	nr := roundNum
	xof := sha3.NewShake256()
	xof.Write(nonce)
	state = make([]uint64, 16)

	rks := make([][]uint64, nr+1)

	for r := 0; r <= nr; r++ {
		rks[r] = make([]uint64, 16)
		for st := 0; st < 16; st++ {
			rks[r][st] = ckks_fv.SampleZqx(xof, plainModulus) * key[st] % plainModulus
		}
	}

	for i := 0; i < 16; i++ {
		state[i] = uint64(i + 1)
	}

	// round0
	for st := 0; st < 16; st++ {
		state[st] = (state[st] + rks[0][st]) % plainModulus
	}

	for r := 1; r < roundNum; r++ {
		for col := 0; col < 4; col++ {
			y0 := 2*state[col] + 3*state[col+4] + 1*state[col+8] + 1*state[col+12]
			y1 := 2*state[col+4] + 3*state[col+8] + 1*state[col+12] + 1*state[col]
			y2 := 2*state[col+8] + 3*state[col+12] + 1*state[col] + 1*state[col+4]
			y3 := 2*state[col+12] + 3*state[col] + 1*state[col+4] + 1*state[col+8]

			state[col] = y0 % plainModulus
			state[col+4] = y1 % plainModulus
			state[col+8] = y2 % plainModulus
			state[col+12] = y3 % plainModulus
		}

		for row := 0; row < 4; row++ {
			y0 := 2*state[4*row] + 3*state[4*row+1] + 1*state[4*row+2] + 1*state[4*row+3]
			y1 := 2*state[4*row+1] + 3*state[4*row+2] + 1*state[4*row+3] + 1*state[4*row]
			y2 := 2*state[4*row+2] + 3*state[4*row+3] + 1*state[4*row] + 1*state[4*row+1]
			y3 := 2*state[4*row+3] + 3*state[4*row] + 1*state[4*row+1] + 1*state[4*row+2]

			state[4*row] = y0 % plainModulus
			state[4*row+1] = y1 % plainModulus
			state[4*row+2] = y2 % plainModulus
			state[4*row+3] = y3 % plainModulus
		}

		for st := 0; st < 16; st++ {
			state[st] = (state[st] * state[st] % plainModulus) * state[st] % plainModulus
		}

		for st := 0; st < 16; st++ {
			state[st] = (state[st] + rks[r][st]) % plainModulus
		}
	}
	for col := 0; col < 4; col++ {
		y0 := 2*state[col] + 3*state[col+4] + 1*state[col+8] + 1*state[col+12]
		y1 := 2*state[col+4] + 3*state[col+8] + 1*state[col+12] + 1*state[col]
		y2 := 2*state[col+8] + 3*state[col+12] + 1*state[col] + 1*state[col+4]
		y3 := 2*state[col+12] + 3*state[col] + 1*state[col+4] + 1*state[col+8]

		state[col] = y0 % plainModulus
		state[col+4] = y1 % plainModulus
		state[col+8] = y2 % plainModulus
		state[col+12] = y3 % plainModulus
	}

	for row := 0; row < 4; row++ {
		y0 := 2*state[4*row] + 3*state[4*row+1] + 1*state[4*row+2] + 1*state[4*row+3]
		y1 := 2*state[4*row+1] + 3*state[4*row+2] + 1*state[4*row+3] + 1*state[4*row]
		y2 := 2*state[4*row+2] + 3*state[4*row+3] + 1*state[4*row] + 1*state[4*row+1]
		y3 := 2*state[4*row+3] + 3*state[4*row] + 1*state[4*row+1] + 1*state[4*row+2]

		state[4*row] = y0 % plainModulus
		state[4*row+1] = y1 % plainModulus
		state[4*row+2] = y2 % plainModulus
		state[4*row+3] = y3 % plainModulus
	}

	for st := 0; st < 16; st++ {
		state[st] = (state[st] * state[st] % plainModulus) * state[st] % plainModulus
	}

	for col := 0; col < 4; col++ {
		y0 := 2*state[col] + 3*state[col+4] + 1*state[col+8] + 1*state[col+12]
		y1 := 2*state[col+4] + 3*state[col+8] + 1*state[col+12] + 1*state[col]
		y2 := 2*state[col+8] + 3*state[col+12] + 1*state[col] + 1*state[col+4]
		y3 := 2*state[col+12] + 3*state[col] + 1*state[col+4] + 1*state[col+8]

		state[col] = y0 % plainModulus
		state[col+4] = y1 % plainModulus
		state[col+8] = y2 % plainModulus
		state[col+12] = y3 % plainModulus
	}

	for row := 0; row < 4; row++ {
		y0 := 2*state[4*row] + 3*state[4*row+1] + 1*state[4*row+2] + 1*state[4*row+3]
		y1 := 2*state[4*row+1] + 3*state[4*row+2] + 1*state[4*row+3] + 1*state[4*row]
		y2 := 2*state[4*row+2] + 3*state[4*row+3] + 1*state[4*row] + 1*state[4*row+1]
		y3 := 2*state[4*row+3] + 3*state[4*row] + 1*state[4*row+1] + 1*state[4*row+2]

		state[4*row] = y0 % plainModulus
		state[4*row+1] = y1 % plainModulus
		state[4*row+2] = y2 % plainModulus
		state[4*row+3] = y3 % plainModulus
	}

	for st := 0; st < 16; st++ {
		state[st] = (state[st] + rks[roundNum][st]) % plainModulus
	}
	return
}
