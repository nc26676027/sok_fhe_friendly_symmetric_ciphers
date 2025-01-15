package main

import (
	"crypto/rand"
	"fmt"
	"math"
	"math/big"
	"time"

	"github.com/ldsec/lattigo/v2/ckks_fv"
	"github.com/ldsec/lattigo/v2/utils"
	"golang.org/x/crypto/sha3"
)

func plainMastaDebug(blocksize int, numRound int, nonce []byte, counter []byte, key []uint64, plainModulus uint64) (state []uint64) {
	xof := sha3.NewShake128()
	xof.Write(nonce)
	xof.Write(counter)
	state = make([]uint64, blocksize)

	rc := make([][]uint64, numRound+1)
	rcMat := make([][][]uint64, numRound+1)

	for r := 0; r <= numRound; r++ {
		rc[r] = make([]uint64, blocksize)
		rc[r] = ckks_fv.GetRandomVector(blocksize, xof, plainModulus)
	}

	for r := 0; r <= numRound; r++ {
		rcMat[r] = make([][]uint64, blocksize)
		rcMat[r] = ckks_fv.GetRandomMatrixMasta(blocksize, xof, plainModulus)
	}

	for i := 0; i < blocksize; i++ {
		key[i] = key[i] % plainModulus
		state[i] = key[i]
	}

	// Round Functions
	for r := 0; r < numRound; r++ {
		mastaLinearLayer(state, rcMat[r], rc[r], plainModulus)
		fmt.Println("After Linlayer state ", r, ", state[10]: ", state[10])

		plainSbox(state, plainModulus)
		fmt.Println("After sbox round: ", r, ", state[10]: ", state[10])
	}
	// Finalization
	mastaLinearLayerDebug(state, rcMat[numRound], rc[numRound], plainModulus)
	fmt.Println("After Linlayer Fin, ", state[10])

	for i := 0; i < blocksize; i++ {
		state[i] = (state[i] + key[i]) % plainModulus
	}
	fmt.Println("After Add Key Fin, ", state[10])
	return
}

func plainMasta(blocksize int, numRound int, nonce []byte, counter []byte, key []uint64, plainModulus uint64) (state []uint64) {
	xof := sha3.NewShake128()
	xof.Write(nonce)
	xof.Write(counter)
	state = make([]uint64, blocksize)

	rc := make([][]uint64, numRound+1)
	rcMat := make([][][]uint64, numRound+1)

	for r := 0; r <= numRound; r++ {
		rc[r] = make([]uint64, blocksize)
		rc[r] = ckks_fv.GetRandomVector(blocksize, xof, plainModulus)
	}
	for r := 0; r <= numRound; r++ {
		rcMat[r] = make([][]uint64, blocksize)
		rcMat[r] = ckks_fv.GetRandomMatrixMasta(blocksize, xof, plainModulus)
	}
	for i := 0; i < blocksize; i++ {
		key[i] = key[i] % plainModulus
		state[i] = key[i]
	}
	// Round Functions
	for r := 0; r < numRound; r++ {
		mastaLinearLayer(state, rcMat[r], rc[r], plainModulus)

		plainSbox(state, plainModulus)
	}
	// Finalization
	mastaLinearLayer(state, rcMat[numRound], rc[numRound], plainModulus)
	for i := 0; i < blocksize; i++ {
		state[i] = (state[i] + key[i]) % plainModulus
	}
	return
}

func mastaLinearLayer(state []uint64, mat [][]uint64, rc []uint64, plainModulus uint64) {
	blocksize := len(state)
	sumPool := make([]uint64, blocksize)
	for row := 0; row < blocksize; row	++{

		buf := make([]uint64, blocksize)
		for col := 0; col < blocksize; col++ {
			tmp := new(big.Int).Mul(new(big.Int).SetUint64(mat[row][col]), new(big.Int).SetUint64(state[col]))
			tmp.Mod(tmp, new(big.Int).SetUint64(plainModulus) )
			buf[col] = tmp.Uint64()
		}
		for i := 1; i < blocksize; i++ {
			buf[0] = (buf[0] + buf[i]) % plainModulus 
		}
		sumPool[row] = buf[0]
	}
	
	for row := 0; row < blocksize; row++ {
		state[row] = (sumPool[row] + rc[row]) % plainModulus
	}

}

func mastaLinearLayerDebug(state []uint64, mat [][]uint64, rc []uint64, plainModulus uint64) {
	blocksize := len(state)
	sumPool := make([]uint64, blocksize)
	for row := 0; row < blocksize; row	++{

		buf := make([]uint64, blocksize)
		for col := 0; col < blocksize; col++ {
			tmp := new(big.Int).Mul(new(big.Int).SetUint64(mat[row][col]), new(big.Int).SetUint64(state[col]))
			tmp.Mod(tmp, new(big.Int).SetUint64(plainModulus) )
			buf[col] = tmp.Uint64()
		}
		for i := 1; i < blocksize; i++ {
			buf[0] = (buf[0] + buf[i]) % plainModulus 
		}
		sumPool[row] = buf[0]
	}
	
	for row := 0; row < blocksize; row++ {
		state[row] = (sumPool[row] + rc[row]) % plainModulus
	}

}


func plainSbox(state []uint64, plainModulus uint64) {
	blocksize := len(state)
	stCp := make([]uint64, 2)
	copy(stCp, state[:2])

	for el := 0; el < blocksize - 2; el++ {
		ind1 := (el + 1) % int(blocksize)
		ind2 := (ind1 + 1) % int(blocksize)

		ind1Big := new(big.Int).SetUint64(state[ind1])
		ind2Big := new(big.Int).SetUint64(state[ind2])
		mult := new(big.Int)
		mult.Mul(ind1Big, ind2Big)
		mult.Mod(mult, new(big.Int).SetUint64(plainModulus) )
		sum := new(big.Int)
		sum.Add(ind2Big, mult)
		sum.Mod(sum, new(big.Int).SetUint64(plainModulus) )

		add := new(big.Int).SetUint64(state[el])
		add.Add(add, sum)
		add.Mod(add, new(big.Int).SetUint64(plainModulus) )
		state[el] = add.Uint64()
	}

	ind1 := blocksize - 1	
	ind2 := 0
	ind1Big := new(big.Int).SetUint64(state[ind1])
	ind2Big := new(big.Int).SetUint64(stCp[ind2])
	mult := new(big.Int)
	mult.Mul(ind1Big, ind2Big)
	mult.Mod(mult, new(big.Int).SetUint64(plainModulus) )
	sum := new(big.Int)
	sum.Add(ind2Big, mult)
	sum.Mod(sum, new(big.Int).SetUint64(plainModulus) )

	add := new(big.Int).SetUint64(state[blocksize - 2])
	add.Add(add, sum)
	add.Mod(add, new(big.Int).SetUint64(plainModulus) )
	state[blocksize - 2] = add.Uint64()

	ind1 = 0
	ind2 = 1
	ind1Big = new(big.Int).SetUint64(stCp[ind1])
	ind2Big = new(big.Int).SetUint64(stCp[ind2])
	mult = new(big.Int)
	mult.Mul(ind1Big, ind2Big)
	mult.Mod(mult, new(big.Int).SetUint64(plainModulus) )
	sum = new(big.Int)
	sum.Add(ind2Big, mult)
	sum.Mod(sum, new(big.Int).SetUint64(plainModulus) )

	add = new(big.Int).SetUint64(state[blocksize - 1])
	add.Add(add, sum)
	add.Mod(add, new(big.Int).SetUint64(plainModulus) )
	state[blocksize - 1] = add.Uint64()

}

func testPlainMasta(mastaParam int) {
	numRound := ckks_fv.MastaParams[mastaParam].NumRound
	blocksize := ckks_fv.MastaParams[mastaParam].Blocksize
	nonce := make([]byte, 8)
	counter := make([]byte, 8)
	key := make([]uint64, blocksize)
	t := ckks_fv. MastaParams[mastaParam].PlainModulus

	// Generate secret key
	for i := 0; i < blocksize; i++ {
		key[i] = uint64(i+1) % t
	}

	rand.Read(nonce)
	rand.Read(counter)
	state := plainMasta(blocksize, numRound, nonce, counter, key, t)
	fmt.Printf("State in Hex: %x\n", state )   //
}

func benchmarkRtFMasta( mastaParam int) {
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
	var masta ckks_fv.MFVMasta

	var data [][]float64
	var nonces [][]byte
	var counter []byte
	var key []uint64
	var keystream [][]uint64
	var fvKeystreams []*ckks_fv.Ciphertext

	// Pasta parameter
	blocksize := ckks_fv.MastaParams[mastaParam].Blocksize
	outputsize := blocksize
	numRound := ckks_fv.MastaParams[mastaParam].NumRound
	plainModulus := ckks_fv.MastaParams[mastaParam].PlainModulus

	// RtF Rubato parameters

	btpParams := ckks_fv.DefaultBootstrapParams[0] // generate btp parameters 
	hbtpParams := btpParams.HalfbtpParams() // get hbtp parameter from btp
	params, err := hbtpParams.Params()
	if err != nil {
		panic(err)
	}
	params.SetPlainModulus(plainModulus)
	params.SetLogFVSlots(params.LogN())
	messageScaling := float64(params.PlainModulus()) / hbtpParams.MessageRatio

	mastaModDown := ckks_fv.MastaModDownParams[mastaParam].CipherModDown
	stcModDown := ckks_fv.MastaModDownParams[mastaParam].StCModDown

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
	_ = plainMastaDebug(blocksize, numRound, nonces[0], counter, key, plainModulus)

	// fmt.Println("Debug plain Done\n\n\n")
	for i := 0; i < params.N(); i++ {
		keystream[i] = plainMasta(blocksize, numRound, nonces[i], counter, key, plainModulus)
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

	masta = ckks_fv.NewMFVMasta( key, mastaParam, params, fvEncoder, fvEncryptor, fvDecryptor, fvEvaluator, mastaModDown[0] )

	Start := time.Now()

	fvKeystreams = masta.Crypt(nonces, counter, mastaModDown)
	for i := 0; i < 1; i++ {
		fmt.Println("Before S2C fvKeystreams level: ", fvKeystreams[i].Level())
		fvKeystreams[i] = fvEvaluator.SlotsToCoeffs(fvKeystreams[i], stcModDown)
		fmt.Println("After S2C fvKeystreams level: ", fvKeystreams[i].Level())
		if fvKeystreams[i].Level() > 0{
			fvEvaluator.ModSwitchMany(fvKeystreams[i], fvKeystreams[i], fvKeystreams[i].Level())
		}
	}
/* We assume that b.N == 1 */
	s2c_time := time.Now()
	ch := make(chan bool, outputsize)
	for i := 0; i < outputsize; i++ {
		go func(i int) {
			evalCopy := fvEvaluator.ShallowCopy()
			fvKeystreams[i] = evalCopy.SlotsToCoeffs(fvKeystreams[i], stcModDown)
			if fvKeystreams[i].Level() > 0{
				fvEvaluator.ModSwitchMany(fvKeystreams[i], fvKeystreams[i], fvKeystreams[i].Level())
			}
			ch <- true
		}(i)
	}
	for i := 0; i < outputsize; i++ {
		<-ch
	}

	// for i := 1; i < outputsize; i++ {
	// 	fvKeystreams[i] = fvEvaluator.SlotsToCoeffs(fvKeystreams[i], stcModDown)
	// 	if fvKeystreams[i].Level() > 0{
	// 		fvEvaluator.ModSwitchMany(fvKeystreams[i], fvKeystreams[i], fvKeystreams[i].Level())
	// 	}
	// }



	s2c_time_dura := time.Since(s2c_time)
	fmt.Printf("The s2c operation took %v\n", s2c_time_dura)

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

func printDebug(params *ckks_fv.Parameters, ciphertext *ckks_fv.Ciphertext, valuesWant []complex128, decryptor ckks_fv.CKKSDecryptor, encoder ckks_fv.CKKSEncoder) {

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
func main() {
	benchmarkRtFMasta( ckks_fv.MASTA7M)
}
