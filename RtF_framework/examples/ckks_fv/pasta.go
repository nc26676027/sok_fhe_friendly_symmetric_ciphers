package main

import (
	"crypto/rand"
	"encoding/binary"
	"fmt"
	"math"
	"math/big"
	"os"
	"time"

	"github.com/ldsec/lattigo/v2/ckks_fv"
	"github.com/ldsec/lattigo/v2/utils"
	"golang.org/x/crypto/sha3"
)



func plainPasta(blocksize int, numRound int, nonce []byte, counter []byte, key []uint64, plainModulus uint64) (state []uint64) {
	xof := sha3.NewShake128()
	xof.Write(nonce)
	xof.Write(counter)
	state = make([]uint64, blocksize)
	PastaT := blocksize / 2

	rc := make([][]uint64, numRound+1)
	rcMat1 := make([][][]uint64, numRound+1)
	rcMat2 := make([][][]uint64, numRound+1)

	for r := 0; r <= numRound; r++ {
		rc[r] = make([]uint64, blocksize)
		rc[r] = ckks_fv.GetRandomVector(blocksize, xof, plainModulus)
	}

	for r := 0; r <= numRound; r++ {
		rcMat1[r] = make([][]uint64, PastaT)
		rcMat2[r] = make([][]uint64, PastaT)
		rcMat1[r] = ckks_fv.GetRandomMatrixPasta(PastaT, xof, plainModulus)
		rcMat2[r] = ckks_fv.GetRandomMatrixPasta(PastaT, xof, plainModulus)
	}

	for i := 0; i < blocksize; i++ {
		state[i] = key[i] % plainModulus
	}

	// // Fisrt Linear layer
	// for i := 0; i < len(state); i++{
	// 	fmt.Println("state[", i, "]", state[i])
	// }
	pastaLinearLayer(state, rcMat1[0], rcMat2[0], rc[0], plainModulus)
	// Round Functions
	for r := 1; r < numRound; r++ {
		pastaFeistel(state, plainModulus)
		pastaLinearLayer(state, rcMat1[r], rcMat2[r], rc[r], plainModulus)
	}
	// Finalization
	for st := 0; st < blocksize; st++ {
		state[st] = (state[st] * state[st] % plainModulus) * state[st] % plainModulus
	}
	// fmt.Println("Plain sboxCube Fin:", state[0])
	pastaLinearLayer(state, rcMat1[numRound], rcMat2[numRound], rc[numRound], plainModulus)
	// fmt.Println("Plain linLayer Fin:", state[0])
	state = state[0 : PastaT]
	return
}

func pastaLinearLayer(state []uint64, mat1, mat2 [][]uint64, rc []uint64, plainModulus uint64) {
	blocksize := len(state)
	PastaT := blocksize / 2
	// fmt.Println("Pasta plain rcmat1: ")
	// ckks_fv.PrintDebugVec(mat1[0], 1, 0)
	// fmt.Println("Pasta plain rcmat2: ")
	// ckks_fv.PrintDebugVec(mat2[0], 1, 0)
	sumPool := make([]uint64, PastaT)
	sumPool2 := make([]uint64, PastaT)
	for row := 0; row < PastaT; row	++{
		buf1 := make([]uint64, PastaT)
		buf2 := make([]uint64, PastaT)

		for col := 0; col < PastaT; col++ {

			tmp := new(big.Int).Mul(new(big.Int).SetUint64(mat1[row][col]), new(big.Int).SetUint64(state[col]))
			tmp.Mod(tmp, new(big.Int).SetUint64(plainModulus) )
			tmp2 := new(big.Int).Mul(new(big.Int).SetUint64(mat2[row][col]), new(big.Int).SetUint64(state[PastaT+col]))
			tmp2.Mod(tmp2, new(big.Int).SetUint64(plainModulus) )
			buf1[col] = tmp.Uint64()
			buf2[col] = tmp2.Uint64()

		}
		// fmt.Println("Plain apply mat1 [", row, "]: ", buf1[10])
		for i := 1; i < PastaT; i++ {
			buf1[0] = (buf1[0] + buf1[i]) % plainModulus 
			buf2[0] = (buf2[0] + buf2[i]) % plainModulus 
		}
		sumPool[row] = buf1[0]
		sumPool2[row] = buf2[0]
	}

	for row := 0; row < PastaT; row++ {
		state[row] = (sumPool[row] + rc[row]) % plainModulus
		state[PastaT+row] = (sumPool2[row] + rc[row+PastaT]) % plainModulus
	}
	for i := 0; i < PastaT; i++{
		sum := ( state[i] + state[PastaT+i] ) % plainModulus
		state[i] = ( state[i] + sum ) % plainModulus
		state[PastaT+i] = ( state[PastaT+i] + sum ) % plainModulus
	}
	// fmt.Println("Pasta plain after mixmat: ")
	// ckks_fv.PrintDebugVec(state, 1, 0)
}

func pastaFeistel(state []uint64, plainModulus uint64) {
	blocksize := len(state)
	for i := blocksize - 1; i > 0; i-- {
		tmp := ( state[i-1] * state[i-1] ) % plainModulus
		state[i] = ( state[i] + tmp ) % plainModulus
	}
	// for i := blocksize - 1; i > 0; i-- {
	// 	fmt.Println("Plain after sboxFeist[", i, "]: ", state[i])
	// }
}

func testPlainPasta(pastaParam int) {
	numRound := ckks_fv.PastaParams[pastaParam].NumRound
	blocksize := ckks_fv.PastaParams[pastaParam].Blocksize
	nonce := make([]byte, 8)
	counter := make([]byte, 8)
	key := make([]uint64, blocksize)
	t := ckks_fv.PastaParams[pastaParam].PlainModulus

	// Generate secret key
	for i := 0; i < blocksize; i++ {
		key[i] = uint64(0) % t
	}

	// Generate nonce
	for i := 0; i < 8; i++ {
		nonce[i] = byte(0)
		counter[i] = byte(0)
	}

	state := plainPasta(blocksize, numRound, nonce, counter, key, t)
	fmt.Printf("State in Hex: %x\n", state )   // 
}

func generateRandomKey(blocksize int) ([]uint64, error) {
	key := make([]uint64, blocksize)
	for i := 0; i < blocksize; i++ {
		var b [8]byte
		if _, err := rand.Read(b[:]); err != nil {
			return nil, err
		}
		key[i] = binary.LittleEndian.Uint64(b[:])
	}
	return key, nil
}

func benchmarkRtFPasta( pastaParam int) {
	var err error

	var hbtp *ckks_fv.HalfBootstrapper
	var kgen ckks_fv.KeyGenerator
	var fvEncoder ckks_fv.MFVEncoder
	var ckksEncoder ckks_fv.CKKSEncoder
	var ckksDecryptor ckks_fv.CKKSDecryptor
	var sk *ckks_fv.SecretKey
	var pk *ckks_fv.PublicKey
	var fvEncryptor ckks_fv.MFVEncryptor
	var fvEvaluator ckks_fv.MFVEvaluator
	var fvDecryptor ckks_fv.MFVDecryptor
	var plainCKKSRingTs []*ckks_fv.PlaintextRingT
	var plaintexts []*ckks_fv.Plaintext
	var pasta ckks_fv.MFVPasta

	var data [][]float64
	var nonces [][]byte
	var counter []byte
	var key []uint64
	var keystream [][]uint64
	var fvKeystreams []*ckks_fv.Ciphertext

	// Pasta parameter
	blocksize := ckks_fv.PastaParams[pastaParam].Blocksize
	outputsize := blocksize / 2
	numRound := ckks_fv.PastaParams[pastaParam].NumRound
	plainModulus := ckks_fv.PastaParams[pastaParam].PlainModulus

	// RtF Rubato parameters
	// var plaintext *ckks_fv.Plaintext
	btpParams := ckks_fv.DefaultBootstrapParams[0] // generate btp parameters 
	hbtpParams := btpParams.HalfbtpParams() // get hbtp parameter from btp

	// hbtpParams := ckks_fv.DefaultBootstrapParams[0]
	params, err := hbtpParams.Params()
	if err != nil {
		panic(err)
	}
	params.SetPlainModulus(plainModulus)
	params.SetLogFVSlots(params.LogN())
	messageScaling := float64(params.PlainModulus()) / hbtpParams.MessageRatio

	pastaModDown := ckks_fv.PastaModDownParams[pastaParam].CipherModDown
	stcModDown := ckks_fv.PastaModDownParams[pastaParam].StCModDown

	// Scheme context and keys
	kgen = ckks_fv.NewKeyGenerator(params)
	sk, pk = kgen.GenKeyPairSparse(hbtpParams.H)
	fvEncoder = ckks_fv.NewMFVEncoder(params)
	ckksEncoder = ckks_fv.NewCKKSEncoder(params)
	fvEncryptor = ckks_fv.NewMFVEncryptorFromPk(params, pk)
	fvDecryptor = ckks_fv.NewMFVDecryptor(params, sk)
	ckksDecryptor = ckks_fv.NewCKKSDecryptor(params, sk)

	// Generating half-bootstrapping keys
	rotationsHalfBoot := kgen.GenRotationIndexesForHalfBoot(params.LogSlots(), hbtpParams)
	pDcds := fvEncoder.GenSlotToCoeffMatFV(2) // radix = 2
	rotationsStC := kgen.GenRotationIndexesForSlotsToCoeffsMat(pDcds)
	rotations := append(rotationsHalfBoot, rotationsStC...)
	rotkeys := kgen.GenRotationKeysForRotations(rotations, true, sk)
	rlk := kgen.GenRelinearizationKey(sk)
	hbtpKey := ckks_fv.BootstrappingKey{Rlk: rlk, Rtks: rotkeys}

	if hbtp, err = ckks_fv.NewHalfBootstrapper(params, hbtpParams, hbtpKey); err != nil {
		panic(err)
	}

	// Encode float data added by keystream to plaintext coefficients
	fvEvaluator = ckks_fv.NewMFVEvaluator(params, ckks_fv.EvaluationKey{Rlk: rlk, Rtks: rotkeys}, pDcds)

	coeffs := make([][]float64, outputsize)
	for s := 0; s < outputsize; s++ {
		coeffs[s] = make([]float64, params.N())
	}

	key = make([]uint64, blocksize)
	for i := 0; i < blocksize; i++ {
		key[i] = uint64(i + 1) % plainModulus // Use (1, ..., 16) for testing
	}

	// Get random data in [-1, 1]
	data = make([][]float64, outputsize)
	for s := 0; s < outputsize; s++ {
		data[s] = make([]float64, params.N())
		for i := 0; i < params.N(); i++ {
			data[s][i] = utils.RandFloat64(-1, 1)
			// data[s][i] = 0
		}
	}

	nonces = make([][]byte, params.N())
	for i := 0; i < params.N(); i++ {
		nonces[i] = make([]byte, 8)
		rand.Read(nonces[i])
	}
	counter = make([]byte, 8)
	rand.Read(counter)

	// Get keystream, focus on the coefficient format of the encoding 
	keystream = make([][]uint64, params.N())
	for i := 0; i < params.N(); i++ {
		keystream[i] = plainPasta(blocksize, numRound, nonces[i], counter, key, plainModulus)
	}

	for s := 0; s < outputsize; s++ {
		for i := 0; i < params.N()/2; i++ {
			j := utils.BitReverse64(uint64(i), uint64(params.LogN()-1))
			coeffs[s][j] = data[s][i]
			coeffs[s][j+uint64(params.N()/2)] = data[s][i+params.N()/2]
		}
	}

	// Encode plaintext
	plainCKKSRingTs = make([]*ckks_fv.PlaintextRingT, outputsize)
	for s := 0; s < outputsize; s++ {
		plainCKKSRingTs[s] = ckksEncoder.EncodeCoeffsRingTNew(coeffs[s], messageScaling)
		poly := plainCKKSRingTs[s].Value()[0]
		for i := 0; i < params.N(); i++ {
			j := utils.BitReverse64(uint64(i), uint64(params.LogN()))
			poly.Coeffs[0][j] = (poly.Coeffs[0][j] + keystream[i][s]) % params.PlainModulus()
		}
	}

	plaintexts = make([]*ckks_fv.Plaintext, outputsize)
	for s := 0; s < outputsize; s++ {
		plaintexts[s] = ckks_fv.NewPlaintextFVLvl(params, 0)
		fvEncoder.FVScaleUp(plainCKKSRingTs[s], plaintexts[s])
	}

	pasta = ckks_fv.NewMFVPasta( key, pastaParam, params, fvEncoder, fvEncryptor, fvDecryptor, fvEvaluator, pastaModDown[0])

	Start := time.Now()

		fvKeystreams = pasta.Crypt(nonces, counter, pastaModDown)
		for i := 0; i < 1; i++ {
			fvKeystreams[i] = fvEvaluator.SlotsToCoeffs(fvKeystreams[i], stcModDown)
			if fvKeystreams[i].Level() > 0 {
				fvEvaluator.ModSwitchMany(fvKeystreams[i], fvKeystreams[i], fvKeystreams[i].Level())
			}
		}
	/* We assume that b.N == 1 */
		for i := 1; i < outputsize; i++ {
			fvKeystreams[i] = fvEvaluator.SlotsToCoeffs(fvKeystreams[i], stcModDown)
			if fvKeystreams[i].Level() > 0{
				fvEvaluator.ModSwitchMany(fvKeystreams[i], fvKeystreams[i], fvKeystreams[i].Level())
			}
			
		}
	elapsed := time.Since(Start) 
	fmt.Printf("The Crypt operation took %v\n", elapsed)

	var ctBoot *ckks_fv.Ciphertext
	// Encrypt and mod switch to the lowest level
	ciphertext := ckks_fv.NewCiphertextFVLvl(params, 1, 0)
	ciphertext.Value()[0] = plaintexts[0].Value()[0].CopyNew()
	fvEvaluator.Sub(ciphertext, fvKeystreams[0], ciphertext)
	fvEvaluator.TransformToNTT(ciphertext, ciphertext)
	ciphertext.SetScale(math.Exp2(math.Round(math.Log2(float64(params.Qi()[0]) / float64(params.PlainModulus()) * messageScaling))))

	// Half-Bootstrap the ciphertext (homomorphic evaluation of ModRaise -> SubSum -> CtS -> EvalMod)
	// It takes a ciphertext at level 0 (if not at level 0, then it will reduce it to level 0)
	// and returns a ciphertext at level MaxLevel - k, where k is the depth of the bootstrapping circuit.
	// Difference from the bootstrapping is that the last StC is missing.
	// CAUTION: the scale of the ciphertext MUST be equal (or very close) to params.Scale
	// To equalize the scale, the function evaluator.SetScale(ciphertext, parameters.Scale) can be used at the expense of one level.
	ctBoot, _ = hbtp.HalfBoot(ciphertext, false)

	elapsed = time.Since(Start) 
	fmt.Printf("The E2E took %v\n", elapsed)
	valuesWant := make([]complex128, params.Slots())

	for i := 0; i < params.Slots(); i++ {
		valuesWant[i] = complex(data[0][i], 0)
	}

	fmt.Println("Precision of HalfBoot(ciphertext)")
	printDebug(params, ctBoot, valuesWant, ckksDecryptor, ckksEncoder)
}


func findPastaModDown(pastaParam int, radix int) {
	var err error

	var kgen ckks_fv.KeyGenerator
	var fvEncoder ckks_fv.MFVEncoder
	var sk *ckks_fv.SecretKey
	var pk *ckks_fv.PublicKey
	var fvEncryptor ckks_fv.MFVEncryptor
	var fvDecryptor ckks_fv.MFVDecryptor
	var fvEvaluator ckks_fv.MFVEvaluator
	var fvNoiseEstimator ckks_fv.MFVNoiseEstimator
	var pasta ckks_fv.MFVPasta

	var nonces [][]byte
	var key []uint64
	var stCt []*ckks_fv.Ciphertext
	var keystream [][]uint64

	var pastaModDown []int
	var stcModDown []int

	// Pasta parameter
	blocksize := ckks_fv.PastaParams[pastaParam].Blocksize
	numRound := ckks_fv.PastaParams[pastaParam].NumRound
	plainModulus := ckks_fv.PastaParams[pastaParam].PlainModulus

	// RtF Pasta parameters
	// Four sets of parameters (index 0 to 1) ensuring 128 bit of security
	// are available in github.com/smilecjf/lattigo/v2/ckks_fv/rtf_params
	// LogSlots is hardcoded in the parameters, but can be changed from 4 to 15.
	// When changing logSlots make sure that the number of levels allocated to CtS is
	// smaller or equal to logSlots.

	btpParams := ckks_fv.DefaultBootstrapParams[0] // generate btp parameters 
	hbtpParams := btpParams.HalfbtpParams() // get hbtp parameter from btp
	params, err := hbtpParams.Params()
	if err != nil {
		panic(err)
	}
	params.SetPlainModulus(plainModulus)
	params.SetLogFVSlots(params.LogN())

	// Scheme context and keys
	kgen = ckks_fv.NewKeyGenerator(params)
	sk, pk = kgen.GenKeyPairSparse(hbtpParams.H)

	fvEncoder = ckks_fv.NewMFVEncoder(params)

	fvEncryptor = ckks_fv.NewMFVEncryptorFromPk(params, pk)
	fvDecryptor = ckks_fv.NewMFVDecryptor(params, sk)
	fvNoiseEstimator = ckks_fv.NewMFVNoiseEstimator(params, sk)

	pDcds := fvEncoder.GenSlotToCoeffMatFV(radix)
	rotations := kgen.GenRotationIndexesForSlotsToCoeffsMat(pDcds)
	rotkeys := kgen.GenRotationKeysForRotations(rotations, true, sk)
	rlk := kgen.GenRelinearizationKey(sk)

	fvEvaluator = ckks_fv.NewMFVEvaluator(params, ckks_fv.EvaluationKey{Rlk: rlk, Rtks: rotkeys}, pDcds)

	// Generating data set
	key = make([]uint64, blocksize)
	// for i := 0; i < blocksize; i++ {
	// 	key[i] = uint64(i + 1) // Use (1, ..., 16) for testing
	// }
	key, _ = generateRandomKey(blocksize)

	nonces = make([][]byte, params.FVSlots())
	for i := 0; i < params.FVSlots(); i++ {
		nonces[i] = make([]byte, 8)
		rand.Read(nonces[i])
	}
	counter := make([]byte, 8)
	rand.Read(counter)

	keystream = make([][]uint64, params.FVSlots())
	for i := 0; i < 1; i++ {
		keystream[i] = plainPasta(blocksize, numRound, nonces[i], counter, key, params.PlainModulus())
	}
	outputsize := blocksize / 2

	Start := time.Now()

	// Find proper nbInitModDown value for fvHera
	fmt.Println("=========== Start to find nbInitModDown ===========")
	pasta = ckks_fv.NewMFVPasta( key, pastaParam, params, fvEncoder, fvEncryptor, fvDecryptor, fvEvaluator, 0)
	stCt = pasta.CryptNoModSwitch(nonces, counter)

	elapsed := time.Since(Start) 
	fmt.Printf("The Crypt operation took %v\n", elapsed)

	invBudgets := make([]int, outputsize)
	minInvBudget := int((^uint(0)) >> 1) // MaxInt
	for i := 0; i < outputsize; i++ {
		ksSlot := fvEvaluator.SlotsToCoeffsNoModSwitch(stCt[i])

		invBudgets[i] = fvNoiseEstimator.InvariantNoiseBudget(ksSlot)
		if invBudgets[i] < minInvBudget {
			minInvBudget = invBudgets[i]
		}
		fvEvaluator.ModSwitchMany(ksSlot, ksSlot, ksSlot.Level())

		ksCt := fvDecryptor.DecryptNew(ksSlot)
		ksCoef := ckks_fv.NewPlaintextRingT(params)
		fvEncoder.DecodeRingT(ksCt, ksCoef)

		for j := 0; j < params.FVSlots(); j++ {
			br_j := utils.BitReverse64(uint64(j), uint64(params.LogN()))

			if ksCoef.Element.Value()[0].Coeffs[0][br_j] != keystream[j][i] {
				fmt.Printf("[-] Validity failed")
				os.Exit(0)
			}
		}
	}
	fmt.Printf("Budget info : min %d in %v\n", minInvBudget, invBudgets)

	qi := params.Qi()
	qiCount := params.QiCount()
	logQi := make([]int, qiCount)
	for i := 0; i < qiCount; i++ {
		logQi[i] = int(math.Round(math.Log2(float64(qi[i]))))
	}

	nbInitModDown := 0
	cutBits := logQi[qiCount-1]
	for cutBits+40 < minInvBudget { // if minInvBudget is too close to cutBits, decryption can be failed
		nbInitModDown++
		cutBits += logQi[qiCount-nbInitModDown-1]
	}
	fmt.Printf("Preferred nbInitModDown = %d\n\n", nbInitModDown)

	fmt.Println("=========== Start to find PastaModDown & StcModDown ===========")
	pasta = ckks_fv.NewMFVPasta(key, pastaParam, params, fvEncoder, fvEncryptor, fvDecryptor, fvEvaluator, nbInitModDown)
	stCt, pastaModDown = pasta.CryptAutoModSwitch(nonces, counter, fvNoiseEstimator)
	_, stcModDown = fvEvaluator.SlotsToCoeffsAutoModSwitch(stCt[0], fvNoiseEstimator)
	for i := 0; i < outputsize; i++ {
		ksSlot := fvEvaluator.SlotsToCoeffs(stCt[i], stcModDown)
		if ksSlot.Level() > 0 {
			fvEvaluator.ModSwitchMany(ksSlot, ksSlot, ksSlot.Level())
		}

		ksCt := fvDecryptor.DecryptNew(ksSlot)
		ksCoef := ckks_fv.NewPlaintextRingT(params)
		fvEncoder.DecodeRingT(ksCt, ksCoef)

		for j := 0; j < params.FVSlots(); j++ {
			br_j := utils.BitReverse64(uint64(j), uint64(params.LogN()))

			if ksCoef.Element.Value()[0].Coeffs[0][br_j] != keystream[j][i] {
				fmt.Printf("[-] Validity failed")
				os.Exit(0)
			}
		}
	}

	fmt.Printf("Pasta modDown : %v\n", pastaModDown)
	fmt.Printf("SlotsToCoeffs modDown : %v\n", stcModDown)
}

func printDebug2(params *ckks_fv.Parameters, ciphertext *ckks_fv.Ciphertext, valuesWant []complex128, decryptor ckks_fv.CKKSDecryptor, encoder ckks_fv.CKKSEncoder) {

	valuesTest := encoder.DecodeComplex(decryptor.DecryptNew(ciphertext), params.LogSlots())
	logSlots := params.LogSlots()
	sigma := params.Sigma()

	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale()))
	fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])

	precStats := ckks_fv.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, logSlots, sigma)

	fmt.Println(precStats.String())
}

func main2() {
	benchmarkRtFPasta( ckks_fv.PASTA5M)
}
