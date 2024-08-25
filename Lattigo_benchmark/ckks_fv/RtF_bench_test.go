package ckks_fv

import (
	"crypto/rand"
	"fmt"
	"math"
	"math/big"
	"testing"
	"time"

	"github.com/ldsec/lattigo/v2/ring"
	"github.com/ldsec/lattigo/v2/utils"
	"golang.org/x/crypto/sha3"
)

// // Benchmark RtF framework with HERA for 80-bit security full-slots parameter
// func BenchmarkRtFHera80f(b *testing.B) {
// 	benchmarkRtFHera(b, "80f", 4, 0, 2, true)
// }

// // Benchmark RtF framework with HERA for 80-bit security 4-slots parameter
// func BenchmarkRtFHera80s(b *testing.B) {
// 	benchmarkRtFHera(b, "80s", 4, 1, 0, false)
// }

// // Benchmark RtF framework with HERA for 80-bit security full-slots parameter with arcsine evaluation
// func BenchmarkRtFHera80af(b *testing.B) {
// 	benchmarkRtFHera(b, "80af", 4, 2, 2, true)
// }

// // Benchmark RtF framework with HERA for 80-bit security 4-slots parameter with arcsine evaluation
// func BenchmarkRtFHera80as(b *testing.B) {
// 	benchmarkRtFHera(b, "80as", 4, 3, 0, false)
// }

// Benchmark RtF framework with HERA for 128-bit security full-slots parameter
func BenchmarkRtFHera128f(b *testing.B) {
	benchmarkRtFHera(b, "128f", 5, 0, 2, true)
}

// // Benchmark RtF framework with HERA for 128-bit security 4-slots parameter
// func BenchmarkRtFHera128s(b *testing.B) {
// 	benchmarkRtFHera(b, "128s", 5, 1, 0, false)
// }

// // Benchmark RtF framework with HERA for 128-bit security full-slots parameter with arcsine evaluation
// func BenchmarkRtFHera128af(b *testing.B) {
// 	benchmarkRtFHera(b, "128af", 5, 2, 2, true)
// }

// // Benchmark RtF framework with HERA for 128-bit security 4-slots parameter with arcsine evaluation
// func BenchmarkRtFHera128as(b *testing.B) {
// 	benchmarkRtFHera(b, "128as", 5, 3, 2, false)
// }

// // Benchmark RtF framework with Rubato80S
// func BenchmarkRtFRubato80S(b *testing.B) {
// 	benchmarkRtFRubato(b, RUBATO80S)
// }

// // Benchmark RtF framework with Rubato80M
// func BenchmarkRtFRubato80M(b *testing.B) {
// 	benchmarkRtFRubato(b, RUBATO80M)
// }

// // Benchmark RtF framework with Rubato80L
// func BenchmarkRtFRubato80L(b *testing.B) {
// 	benchmarkRtFRubato(b, RUBATO80L)
// }

// Benchmark RtF framework with Rubato128S
func BenchmarkRtFRubato128S(b *testing.B) {
	benchmarkRtFRubato(b, RUBATO128S)
}

// Benchmark RtF framework with Rubato128M
func BenchmarkRtFRubato128M(b *testing.B) {
	benchmarkRtFRubato(b, RUBATO128M)
}

// Benchmark RtF framework with Rubato128L
func BenchmarkRtFRubato128L(b *testing.B) {
	benchmarkRtFRubato(b, RUBATO128L)
}

// // Benchmark RtF framework with HERA for 128-bit security 4-slots parameter
// func BenchmarkRtFPasta4M(b *testing.B) {
// 	benchmarkRtFPasta(b, PASTA4M)
// }


// benchmarkRtFHera(b, "80f", 4, 0, 2, true)
func benchmarkRtFHera(b *testing.B, name string, numRound int, paramIndex int, radix int, fullCoeffs bool) {
	var err error

	var hbtp *HalfBootstrapper
	var kgen KeyGenerator
	var fvEncoder MFVEncoder
	var ckksEncoder CKKSEncoder
	var ckksDecryptor CKKSDecryptor
	var sk *SecretKey
	var pk *PublicKey
	var fvEncryptor MFVEncryptor
	var fvEvaluator MFVEvaluator
	var plainCKKSRingTs []*PlaintextRingT
	var plaintexts []*Plaintext
	var hera MFVHera

	var data [][]float64
	var nonces [][]byte
	var key []uint64
	var keystream [][]uint64
	var fvKeystreams []*Ciphertext

	// RtF Hera parameters
	// Four sets of parameters (index 0 to 3) ensuring 128 bit of security
	// are available in github.com/smilecjf/lattigo/v2/ckks_fv/rtf_params
	// LogSlots is hardcoded in the parameters, but can be changed from 4 to 15.
	// When changing logSlots make sure that the number of levels allocated to CtS is
	// smaller or equal to logSlots.

	// hbtpParams := RtFHeraParams[paramIndex]
	btpParams := DefaultBootstrapParams[0]
	hbtpParams := btpParams.HalfbtpParams()
	params, err := hbtpParams.Params()
	if err != nil {
		panic(err)
	}
	messageScaling := float64(params.PlainModulus()) / hbtpParams.MessageRatio

	// HERA parameters in RtF
	var heraModDown, stcModDown []int
	if numRound == 4 {
		heraModDown = HeraModDownParams80[paramIndex].CipherModDown
		stcModDown = HeraModDownParams80[paramIndex].StCModDown
	} else {
		heraModDown = HeraModDownParams128[paramIndex].CipherModDown
		stcModDown = HeraModDownParams128[paramIndex].StCModDown
	}

	// fullCoeffs denotes whether full coefficients are used for data encoding
	if fullCoeffs {
		params.SetLogFVSlots(params.LogN())
	} else {
		params.SetLogFVSlots(params.LogSlots())
	}

	// Scheme context and keys
	kgen = NewKeyGenerator(params)

	sk, pk = kgen.GenKeyPairSparse(hbtpParams.H)

	fvEncoder = NewMFVEncoder(params)
	ckksEncoder = NewCKKSEncoder(params)
	fvEncryptor = NewMFVEncryptorFromPk(params, pk)
	ckksDecryptor = NewCKKSDecryptor(params, sk)

	// Generating half-bootstrapping keys
	rotationsHalfBoot := kgen.GenRotationIndexesForHalfBoot(params.LogSlots(), hbtpParams)
	pDcds := fvEncoder.GenSlotToCoeffMatFV(radix)
	rotationsStC := kgen.GenRotationIndexesForSlotsToCoeffsMat(pDcds)
	rotations := append(rotationsHalfBoot, rotationsStC...)
	if !fullCoeffs {
		rotations = append(rotations, params.Slots()/2)
	}
	rotkeys := kgen.GenRotationKeysForRotations(rotations, true, sk)
	rlk := kgen.GenRelinearizationKey(sk)
	hbtpKey := BootstrappingKey{Rlk: rlk, Rtks: rotkeys}

	if hbtp, err = NewHalfBootstrapper(params, hbtpParams, hbtpKey); err != nil {
		panic(err)
	}

	// Encode float data added by keystream to plaintext coefficients
	fvEvaluator = NewMFVEvaluator(params, EvaluationKey{Rlk: rlk, Rtks: rotkeys}, pDcds)
	
	coeffs := make([][]float64, 16)
	for s := 0; s < 16; s++ { // 16 polynomials to store 16block of Hera, using row packing 
		coeffs[s] = make([]float64, params.N())
	}

	key = make([]uint64, 16)
	for i := 0; i < 16; i++ {
		key[i] = uint64(i + 1) // Use (1, ..., 16) for testing
	}

	if fullCoeffs {
		fmt.Println("fullSlots, function")
		data = make([][]float64, 16)
		for s := 0; s < 16; s++ {
			data[s] = make([]float64, params.N())
			for i := 0; i < params.N() / 2; i++ {
				data[s][i] = utils.RandFloat64(-1, 1)
				// data[s][i] = 0.5
			}
		}

		nonces = make([][]byte, params.N())
		for i := 0; i < params.N(); i++ {
			nonces[i] = make([]byte, 64)
			rand.Read(nonces[i])
		}

		keystream = make([][]uint64, params.N())
		for i := 0; i < params.N(); i++ {
			keystream[i] = plainHera(numRound, nonces[i], key, params.PlainModulus())
		}
		// fmt.Println("total parmas.N block of symmetric enc, keystream block 0: ", keystream[0])

		for s := 0; s < 16; s++ {
			for i := 0; i < params.N()/2; i++ {
				j := utils.BitReverse64(uint64(i), uint64(params.LogN()-1))
				coeffs[s][j] = data[s][i] // ?? first half data encoded to front part coeff
				coeffs[s][j+uint64(params.N()/2)] = data[s][i+params.N()/2]// ?? second half data encoded to behand part coeff 
			}
		}

		plainCKKSRingTs = make([]*PlaintextRingT, 16)
		for s := 0; s < 16; s++ {
			plainCKKSRingTs[s] = ckksEncoder.EncodeCoeffsRingTNew(coeffs[s], messageScaling) //encode plain into coeff
			poly := plainCKKSRingTs[s].Value()[0] // RNS of RingTs has only one moduli
			for i := 0; i < params.N(); i++ {
				j := utils.BitReverse64(uint64(i), uint64(params.LogN()))
				poly.Coeffs[0][j] = (poly.Coeffs[0][j] + keystream[i][s]) % params.PlainModulus()
			}
		}
	} else {
		fmt.Println("no fullSlots: ")
		data = make([][]float64, 16)
		for s := 0; s < 16; s++ {
			data[s] = make([]float64, params.Slots())
			for i := 0; i < params.Slots(); i++ {
				data[s][i] = utils.RandFloat64(-1, 1)// Read random data 
			}
		}

		nonces = make([][]byte, params.Slots())
		for i := 0; i < params.Slots(); i++ {
			nonces[i] = make([]byte, 64)
			rand.Read(nonces[i])
		}

		keystream = make([][]uint64, params.Slots())
		for i := 0; i < params.Slots(); i++ {
			keystream[i] = plainHera(numRound, nonces[i], key, params.PlainModulus())
		}

		for s := 0; s < 16; s++ {
			for i := 0; i < params.Slots()/2; i++ {
				j := utils.BitReverse64(uint64(i), uint64(params.LogN()-1))
				coeffs[s][j] = data[s][i]
				coeffs[s][j+uint64(params.N()/2)] = data[s][i+params.Slots()/2]
			}
		}

		plainCKKSRingTs = make([]*PlaintextRingT, 16)
		for s := 0; s < 16; s++ {
			plainCKKSRingTs[s] = ckksEncoder.EncodeCoeffsRingTNew(coeffs[s], messageScaling)
			poly := plainCKKSRingTs[s].Value()[0]
			for i := 0; i < params.Slots(); i++ {
				j := utils.BitReverse64(uint64(i), uint64(params.LogN()))
				poly.Coeffs[0][j] = (poly.Coeffs[0][j] + keystream[i][s]) % params.PlainModulus()
			}
		}

	}

	plaintexts = make([]*Plaintext, 16) // plaintext used to store plainCKKSRingTs obtained from symmetric

	for s := 0; s < 16; s++ {
		plaintexts[s] = NewPlaintextFVLvl(params, 0) // generate a level one FV plaintext
		fvEncoder.FVScaleUp(plainCKKSRingTs[s], plaintexts[s]) // transforms a PlaintextRingT (R_t) into a Plaintext (R_q) by scaling up the coefficient by Q/t.
	}

	hera = NewMFVHera(numRound, params, fvEncoder, fvEncryptor, fvEvaluator, heraModDown[0]) // initial hera decryption object
	kCt := hera.EncKey(key) // obtion the FV encrypted symmetric key of hera
	Start := time.Now()
	// FV Keystream
	benchOffLat := fmt.Sprintf("RtF HERA Offline Latency")
	b.Run(benchOffLat, func(b *testing.B) {
		fvKeystreams = hera.Crypt(nonces, kCt, heraModDown) // homomorphically generate keystream encrypted under FV
		for i := 0; i < 1; i++ {
			fvKeystreams[i] = fvEvaluator.SlotsToCoeffs(fvKeystreams[i], stcModDown) // FV S2C to get Coeff format keystream
			if fvKeystreams[i].Level() > 0{
				fvEvaluator.ModSwitchMany(fvKeystreams[i], fvKeystreams[i], fvKeystreams[i].Level())
			}
		}
	})
	elapsedTime := time.Since(Start)
	fmt.Print("Time Crypt and Halfboot first ct: ", elapsedTime)

	/* We assume that b.N == 1 */
	benchOffThrput := fmt.Sprintf("RtF HERA Offline Throughput")
	b.Run(benchOffThrput, func(b *testing.B) {
		for i := 1; i < 16; i++ {
			fvKeystreams[i] = fvEvaluator.SlotsToCoeffs(fvKeystreams[i], stcModDown) // homomorphically generate keystream encrypted under FV
			if fvKeystreams[i].Level() > 0{
				fvEvaluator.ModSwitchMany(fvKeystreams[i], fvKeystreams[i], fvKeystreams[i].Level())
			}
		}
	})


	var ctBoot, ct2 *Ciphertext
	benchOnline := fmt.Sprintf("RtF HERA Online Lat x1")
	b.Run(benchOnline, func(b *testing.B) {
		// Encrypt and mod switch to the lowest level
		ciphertext := NewCiphertextFVLvl(params, 1, 0)
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
		if fullCoeffs {
			ctBoot, ct2 = hbtp.HalfBoot(ciphertext, false)
		} else {
			ctBoot, ct2 = hbtp.HalfBoot(ciphertext, true)
		}
	})
	valuesWant := make([]complex128, params.Slots())
	for i := 0; i < params.Slots(); i++ {
		valuesWant[i] = complex(data[0][i], 0)
	}
	valuesWant2 := make([]complex128, params.Slots())
	for i := params.Slots(); i < params.N(); i++ {
		valuesWant2[i-params.Slots()] = complex(data[0][i], 0)
	}
	fmt.Println("before operate ct0")
	printDebug(params, ctBoot, valuesWant, ckksDecryptor, ckksEncoder)
	fmt.Println("before operate ct2")
	printDebug(params, ct2, valuesWant2, ckksDecryptor, ckksEncoder)

	hbtp.ckksEvaluator.Mul(ctBoot, ctBoot, ctBoot)
	hbtp.ckksEvaluator.Rescale(ctBoot, float64(ctBoot.Level()), ctBoot)
	fmt.Println("Precision of HalfBoot(ciphertext)")
	printDebug(params, ctBoot, valuesWant, ckksDecryptor, ckksEncoder)
}

func benchmarkRtFRubato(b *testing.B, rubatoParam int) {
	var err error

	var hbtp *HalfBootstrapper
	var kgen KeyGenerator
	var fvEncoder MFVEncoder
	var ckksEncoder CKKSEncoder
	var ckksDecryptor CKKSDecryptor
	var sk *SecretKey
	var pk *PublicKey
	var fvEncryptor MFVEncryptor
	var fvEvaluator MFVEvaluator
	var plainCKKSRingTs []*PlaintextRingT
	var plaintexts []*Plaintext
	var rubato MFVRubato

	var data [][]float64
	var nonces [][]byte
	var counter []byte
	var key []uint64
	var keystream [][]uint64
	var fvKeystreams []*Ciphertext

	// Rubato parameter
	blocksize := RubatoParams[rubatoParam].Blocksize
	outputsize := blocksize - 4
	numRound := RubatoParams[rubatoParam].NumRound
	plainModulus := RubatoParams[rubatoParam].PlainModulus
	sigma := RubatoParams[rubatoParam].Sigma

	// RtF Rubato parameters
	// hbtpParams := RtFRubatoParams[0]
	btpParams := DefaultBootstrapParams[0]
	hbtpParams := btpParams.HalfbtpParams()
	params, err := hbtpParams.Params()
	if err != nil {
		panic(err)
	}
	params.SetPlainModulus(plainModulus)
	params.SetLogFVSlots(params.LogN())
	messageScaling := float64(params.PlainModulus()) / hbtpParams.MessageRatio

	rubatoModDown := RubatoModDownParams[rubatoParam].CipherModDown
	stcModDown := RubatoModDownParams[rubatoParam].StCModDown

	// Scheme context and keys
	kgen = NewKeyGenerator(params)
	sk, pk = kgen.GenKeyPairSparse(hbtpParams.H)

	fvEncoder = NewMFVEncoder(params)
	ckksEncoder = NewCKKSEncoder(params)
	fvEncryptor = NewMFVEncryptorFromPk(params, pk)
	ckksDecryptor = NewCKKSDecryptor(params, sk)

	// Generating half-bootstrapping keys
	rotationsHalfBoot := kgen.GenRotationIndexesForHalfBoot(params.LogSlots(), hbtpParams)
	pDcds := fvEncoder.GenSlotToCoeffMatFV(2) // radix = 2
	rotationsStC := kgen.GenRotationIndexesForSlotsToCoeffsMat(pDcds)
	rotations := append(rotationsHalfBoot, rotationsStC...)
	rotkeys := kgen.GenRotationKeysForRotations(rotations, true, sk)
	rlk := kgen.GenRelinearizationKey(sk)
	hbtpKey := BootstrappingKey{Rlk: rlk, Rtks: rotkeys}

	if hbtp, err = NewHalfBootstrapper(params, hbtpParams, hbtpKey); err != nil {
		panic(err)
	}

	// Encode float data added by keystream to plaintext coefficients
	fvEvaluator = NewMFVEvaluator(params, EvaluationKey{Rlk: rlk, Rtks: rotkeys}, pDcds)
	coeffs := make([][]float64, outputsize)
	for s := 0; s < outputsize; s++ {
		coeffs[s] = make([]float64, params.N())
	}

	key = make([]uint64, blocksize)
	for i := 0; i < blocksize; i++ {
		key[i] = uint64(i + 1) // Use (1, ..., 16) for testing
	}

	// Get random data in [-1, 1]
	data = make([][]float64, outputsize)
	for s := 0; s < outputsize; s++ {
		data[s] = make([]float64, params.N())
		for i := 0; i < params.N(); i++ {
			data[s][i] = utils.RandFloat64(-1, 1)
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
		keystream[i] = plainRubato(blocksize, numRound, nonces[i], counter, key, plainModulus, sigma)
	}

	for s := 0; s < outputsize; s++ {
		for i := 0; i < params.N()/2; i++ {
			j := utils.BitReverse64(uint64(i), uint64(params.LogN()-1))
			coeffs[s][j] = data[s][i]
			coeffs[s][j+uint64(params.N()/2)] = data[s][i+params.N()/2]
		}
	}

	// Encode plaintext
	plainCKKSRingTs = make([]*PlaintextRingT, outputsize)
	for s := 0; s < outputsize; s++ {
		plainCKKSRingTs[s] = ckksEncoder.EncodeCoeffsRingTNew(coeffs[s], messageScaling)
		poly := plainCKKSRingTs[s].Value()[0]
		for i := 0; i < params.N(); i++ {
			j := utils.BitReverse64(uint64(i), uint64(params.LogN()))
			poly.Coeffs[0][j] = (poly.Coeffs[0][j] + keystream[i][s]) % params.PlainModulus()
		}
	}

	plaintexts = make([]*Plaintext, outputsize)
	for s := 0; s < outputsize; s++ {
		plaintexts[s] = NewPlaintextFVLvl(params, 0)
		fvEncoder.FVScaleUp(plainCKKSRingTs[s], plaintexts[s])
	}

	// FV Keystream
	rubato = NewMFVRubato(rubatoParam, params, fvEncoder, fvEncryptor, fvEvaluator, rubatoModDown[0])
	kCt := rubato.EncKey(key)

	Start := time.Now()
	benchOffLat := fmt.Sprintf("RtF Rubato Offline Latency")
	b.Run(benchOffLat, func(b *testing.B) {
		fvKeystreams = rubato.Crypt(nonces, counter, kCt, rubatoModDown)
		for i := 0; i < 1; i++ {
			fvKeystreams[i] = fvEvaluator.SlotsToCoeffs(fvKeystreams[i], stcModDown)
			fvEvaluator.ModSwitchMany(fvKeystreams[i], fvKeystreams[i], fvKeystreams[i].Level())
		}
	})
	elapsedTime := time.Since(Start)
	fmt.Print("Time Crypt and Halfboot first ct: ", elapsedTime)
	/* We assume that b.N == 1 */
	benchOffThrput := fmt.Sprintf("RtF Rubato Offline Throughput")
	b.Run(benchOffThrput, func(b *testing.B) {
		for i := 1; i < outputsize; i++ {
			fvKeystreams[i] = fvEvaluator.SlotsToCoeffs(fvKeystreams[i], stcModDown)
			fvEvaluator.ModSwitchMany(fvKeystreams[i], fvKeystreams[i], fvKeystreams[i].Level())
		}
	})

	var ctBoot *Ciphertext
	benchOnline := fmt.Sprintf("RtF Rubato Online Lat x1")
	b.Run(benchOnline, func(b *testing.B) {
		// Encrypt and mod switch to the lowest level
		ciphertext := NewCiphertextFVLvl(params, 1, 0)
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
	})
	valuesWant := make([]complex128, params.Slots())
	for i := 0; i < params.Slots(); i++ {
		valuesWant[i] = complex(data[0][i], 0)
	}

	fmt.Println("Precision of HalfBoot(ciphertext)")
	printDebug(params, ctBoot, valuesWant, ckksDecryptor, ckksEncoder)
}


func benchmarkRtFPasta(b *testing.B, pastaParam int) {
	var err error

	var hbtp *HalfBootstrapper
	var kgen KeyGenerator
	var fvEncoder MFVEncoder
	var ckksEncoder CKKSEncoder
	var ckksDecryptor CKKSDecryptor
	var sk *SecretKey
	var pk *PublicKey
	var fvEncryptor MFVEncryptor
	var fvEvaluator MFVEvaluator
	var plainCKKSRingTs []*PlaintextRingT
	var plaintexts []*Plaintext
	var pasta MFVPasta

	var data [][]float64
	var nonces [][]byte
	var counter []byte
	var key []uint64
	var keystream [][]uint64
	var fvKeystreams []*Ciphertext

	// Pasta parameter
	blocksize := PastaParams[pastaParam].Blocksize
	outputsize := blocksize / 2
	numRound := PastaParams[pastaParam].NumRound
	plainModulus := PastaParams[pastaParam].PlainModulus

	// RtF Rubato parameters
	// hbtpParams := RtFHeraParams[0]
	btpParams := DefaultBootstrapParams[0]
	hbtpParams := btpParams.HalfbtpParams()
	params, err := hbtpParams.Params()
	if err != nil {
		panic(err)
	}
	params.SetPlainModulus(plainModulus)
	params.SetLogFVSlots(params.LogN())
	messageScaling := float64(params.PlainModulus()) / hbtpParams.MessageRatio

	pastaModDown := PastaModDownParams[pastaParam].CipherModDown
	stcModDown := PastaModDownParams[pastaParam].StCModDown

	// Scheme context and keys
	kgen = NewKeyGenerator(params)
	sk, pk = kgen.GenKeyPairSparse(hbtpParams.H)
	fvEncoder = NewMFVEncoder(params)
	ckksEncoder = NewCKKSEncoder(params)
	fvEncryptor = NewMFVEncryptorFromPk(params, pk)
	fvDecryptor := NewMFVDecryptor(params, sk)
	ckksDecryptor = NewCKKSDecryptor(params, sk)

	// Generating half-bootstrapping keys
	rotationsHalfBoot := kgen.GenRotationIndexesForHalfBoot(params.LogSlots(), hbtpParams)
	pDcds := fvEncoder.GenSlotToCoeffMatFV(2) // radix = 2
	rotationsStC := kgen.GenRotationIndexesForSlotsToCoeffsMat(pDcds)
	rotations := append(rotationsHalfBoot, rotationsStC...)
	rotkeys := kgen.GenRotationKeysForRotations(rotations, true, sk)
	rlk := kgen.GenRelinearizationKey(sk)
	hbtpKey := BootstrappingKey{Rlk: rlk, Rtks: rotkeys}

	if hbtp, err = NewHalfBootstrapper(params, hbtpParams, hbtpKey); err != nil {
		panic(err)
	}

	// Encode float data added by keystream to plaintext coefficients
	fvEvaluator = NewMFVEvaluator(params, EvaluationKey{Rlk: rlk, Rtks: rotkeys}, pDcds)

	coeffs := make([][]float64, outputsize)
	for s := 0; s < outputsize; s++ {
		coeffs[s] = make([]float64, params.N())
	}

	key = make([]uint64, blocksize)
	for i := 0; i < blocksize; i++ {
		key[i] = uint64(i + 1) // Use (1, ..., 16) for testing
	}

	// Get random data in [-1, 1]
	data = make([][]float64, outputsize)
	for s := 0; s < outputsize; s++ {
		data[s] = make([]float64, params.N())
		for i := 0; i < params.N(); i++ {
			data[s][i] = utils.RandFloat64(-1, 1)
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
	plainCKKSRingTs = make([]*PlaintextRingT, outputsize)
	for s := 0; s < outputsize; s++ {
		plainCKKSRingTs[s] = ckksEncoder.EncodeCoeffsRingTNew(coeffs[s], messageScaling)
		poly := plainCKKSRingTs[s].Value()[0]
		for i := 0; i < params.N(); i++ {
			j := utils.BitReverse64(uint64(i), uint64(params.LogN()))
			poly.Coeffs[0][j] = (poly.Coeffs[0][j] + keystream[i][s]) % params.PlainModulus()
		}
	}

	plaintexts = make([]*Plaintext, outputsize)
	for s := 0; s < outputsize; s++ {
		plaintexts[s] = NewPlaintextFVLvl(params, 0)
		fvEncoder.FVScaleUp(plainCKKSRingTs[s], plaintexts[s])
	}


	pasta = NewMFVPasta( key, pastaParam, params, fvEncoder, fvEncryptor, fvDecryptor, fvEvaluator, pastaModDown[0])

	benchOffLat := fmt.Sprintf("RtF Pasta Offline Latency")
	b.Run(benchOffLat, func(b *testing.B) {
		fvKeystreams = pasta.Crypt(nonces, counter, pastaModDown)
		for i := 0; i < 1; i++ {
			fvKeystreams[i] = fvEvaluator.SlotsToCoeffs(fvKeystreams[i], stcModDown)
			fvEvaluator.ModSwitchMany(fvKeystreams[i], fvKeystreams[i], fvKeystreams[i].Level())
		}
	})
	/* We assume that b.N == 1 */
	benchOffThrput := fmt.Sprintf("RtF Pasta Offline Throughput")
	b.Run(benchOffThrput, func(b *testing.B) {
		for i := 1; i < outputsize; i++ {
			fvKeystreams[i] = fvEvaluator.SlotsToCoeffs(fvKeystreams[i], stcModDown)
			fvEvaluator.ModSwitchMany(fvKeystreams[i], fvKeystreams[i], fvKeystreams[i].Level())
		}
	})

	var ctBoot *Ciphertext
	benchOnline := fmt.Sprintf("RtF Pasta Online Lat x1")
	b.Run(benchOnline, func(b *testing.B) {
		// Encrypt and mod switch to the lowest level
		ciphertext := NewCiphertextFVLvl(params, 1, 0)
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
	})
	valuesWant := make([]complex128, params.Slots())
	for i := 0; i < params.Slots(); i++ {
		valuesWant[i] = complex(data[0][i], 0)
	}

	fmt.Println("Precision of HalfBoot(ciphertext)")
	printDebug(params, ctBoot, valuesWant, ckksDecryptor, ckksEncoder)
}

func printDebug(params *Parameters, ciphertext *Ciphertext, valuesWant []complex128, decryptor CKKSDecryptor, encoder CKKSEncoder) {

	valuesTest := encoder.DecodeComplex(decryptor.DecryptNew(ciphertext), params.LogSlots())
	logSlots := params.LogSlots()
	sigma := params.Sigma()

	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))
	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale()))
	fmt.Printf("ValuesTest: %6.10f %6.10f %6.10f %6.10f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %6.10f %6.10f %6.10f %6.10f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])

	precStats := GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, logSlots, sigma)

	fmt.Println(precStats.String())
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
			rks[r][st] = SampleZqx(xof, plainModulus) * key[st] % plainModulus
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
			rks[r][i] = SampleZqx(xof, plainModulus) * key[i] % plainModulus
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
		rc[r] = GetRandomVector(blocksize, xof, plainModulus)
	}

	for r := 0; r <= numRound; r++ {
		rcMat1[r] = make([][]uint64, PastaT)
		rcMat2[r] = make([][]uint64, PastaT)
		rcMat1[r] = GetRandomMatrixPasta(PastaT, xof, plainModulus)
		rcMat2[r] = GetRandomMatrixPasta(PastaT, xof, plainModulus)
	}

	for i := 0; i < blocksize; i++ {
		state[i] = key[i]
	}

	// Fisrt Linear layer
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
	pastaLinearLayer(state, rcMat1[numRound], rcMat2[numRound], rc[numRound], plainModulus)

	state = state[0 : PastaT]
	return
}

func pastaLinearLayer(state []uint64, mat1, mat2 [][]uint64, rc []uint64, plainModulus uint64) {
	blocksize := len(state)
	PastaT := blocksize / 2
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
		for i := 1; i < len(buf1); i++ {
			buf1[0] = buf1[0] + buf1[i]
			buf2[0] = buf2[0] + buf2[i]
		}
		buf1[0] = buf1[0] + rc[row]
		buf2[0] = buf2[0] + rc[row]
		state[row] = buf1[0] % plainModulus
		state[PastaT+row] = buf2[0] % plainModulus
	}
	for i := 0; i < PastaT; i++{
		sum := ( state[i] + state[PastaT+i] ) % plainModulus
		state[i] = ( state[i] + sum ) % plainModulus
		state[PastaT+i] = ( state[PastaT+i] + sum ) % plainModulus
	}

}

func pastaFeistel(state []uint64, plainModulus uint64) {
	blocksize := len(state)
	buf := make([]uint64, blocksize)

	for i := 0; i < blocksize; i++ {
		buf[i] = state[i]
	}

	for i := 1; i < blocksize; i++ {
		state[i] = (buf[i] + buf[i-1]*buf[i-1]) % plainModulus
	}
}


func testPlainPasta(pastaParam int) {
	numRound := PastaParams[pastaParam].NumRound
	blocksize := PastaParams[pastaParam].Blocksize
	nonce := make([]byte, 8)
	counter := make([]byte, 8)
	key := make([]uint64, blocksize)
	t := PastaParams[pastaParam].PlainModulus

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
	fmt.Printf("State in Hex: %x\n", state )   // 小写state)
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
