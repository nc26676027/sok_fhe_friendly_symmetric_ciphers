package minimax_comp

import (
	"crypto/rand"
	"fmt"
	"math"
	"time"

	"github.com/ldsec/lattigo/v2/ckks_fv"
	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"golang.org/x/crypto/sha3"
)

func MinimaxCompRubato( name string, rubatoParam int) {


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

	galStepsVector := []int{0}
	for i:=0;i<context.params.LogN()-1;i++{
		galStepsVector = append(galStepsVector, 1<<i)
	}
	rotations := []int { }
	galStepsVector = UnitSlice(galStepsVector, rotations)
	galStepsVector = UnitSlice(galStepsVector, rotationsRtF) //get RtF rotation key
	// sort.Ints(rotations) // rotations indices must under order
	fmt.Println("Repeated free rotations: ", galStepsVector)
	rotkeys := context.kgen.GenRotationKeysForRotations(galStepsVector, true, context.sk)

	evk := ckks_fv.EvaluationKey{Rlk: context.rlk, Rtks: rotkeys}

	fvEvaluator = ckks_fv.NewMFVEvaluator(context.params, evk, pDcds)
	context.evaluator = ckks_fv.NewCKKSEvaluator(context.params, evk)
	var hbtp *ckks_fv.HalfBootstrapper

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

	for s := 0; s < outputsize; s++ {
		data[s] = make([]float64, context.params.N())
		for i:=0;i<len(data[s]);i++ {
			if i< (context.params.N()/2){
				data[s][i] = 0.45
			} else{
				data[s][i] = 0.8
			}
		}
	}

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

	// panic("Done RtF framework")
	ResidualModuli := ckks_fv.DefaultBootstrapParams[0].ResidualModuli
	fmt.Println("totalLevel", context.params.QiCount(), "ResidualModuli: ",len( ResidualModuli) )
	context.remainLevel = len( ResidualModuli ) - 1
	logn := context.params.LogSlots()
	fmt.Println("LogSlots: ", logn)

	Start := time.Now()
	ctxt, _ := minimaxComp(comp_no, deg, alpha, tree, scaledVal, context, ctBoot, ctBoot2)

	rtnVec := context.encoder.DecodeComplex(context.decryptor.DecryptNew(ctxt), context.params.LogSlots())
	ckks_fv.PrintDebugVecC(rtnVec, 7, 3)
	elapsed := time.Since(Start) // 
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
