package main

import (
	"flag"
	"fmt"
	"runtime"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/ckks_cipher"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

var flagShort = flag.Bool("short", false, "run the example with a smaller and insecure ring degree.")


func main() {

	flag.Parse()
	runtime.GOMAXPROCS(64)
	// Default LogN, which with the following defined parameters
	// provides a security of 128-bit.

	LogN := 15
	LogDefaultScale := 35

	q0 := []int{36}                                    // 3) ScaleDown & 4) ModUp
	qiSlotsToCoeffs := []int{35, 35}               // 1) SlotsToCoeffs
	qiCircuitSlots := []int{35, 35, 35}           // 0) Circuit in the slot domain
	qiEvalMod := []int{36, 36, 36, 36, 36, 36, 36, 36} // 6) EvalMod
	qiCoeffsToSlots := []int{36, 36, 36}           // 5) CoeffsToSlots

	LogQ := append(q0, qiSlotsToCoeffs...)
	LogQ = append(LogQ, qiCircuitSlots...)
	LogQ = append(LogQ, qiEvalMod...)
	LogQ = append(LogQ, qiCoeffsToSlots...)

	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            LogN,                                              // Log2 of the ring degree
		LogQ:            LogQ, // Log2 of the ciphertext prime moduli
		LogP:            []int{36, 36, 36, 36, 36},                                 // Log2 of the key-switch auxiliary prime moduli
		LogDefaultScale: LogDefaultScale,                                                // Log2 of the scale
		Xs:              ring.Ternary{H: 192},
	})
	if err != nil {
		panic(err)
	}

	// CoeffsToSlots parameters (homomorphic encoding)
	CoeffsToSlotsParameters := dft.MatrixLiteral{
		Type:         dft.HomomorphicEncode,
		Format:       dft.RepackImagAsReal, // Returns the real and imaginary part into separate ciphertexts
		LogSlots:     params.LogMaxSlots(),
		LevelQ:       params.MaxLevelQ(),
		LevelP:       params.MaxLevelP(),
		LogBSGSRatio: 1,
		Levels:       []int{1, 1, 1}, //qiCoeffsToSlots
	}

	// Parameters of the homomorphic modular reduction x mod 1
	Mod1ParametersLiteral := mod1.ParametersLiteral{
		LevelQ:          params.MaxLevel() - CoeffsToSlotsParameters.Depth(true),
		LogScale:        36,               // Matches qiEvalMod
		Mod1Type:        mod1.CosDiscrete, // Multi-interval Chebyshev interpolation
		Mod1Degree:      254,               // Depth 5
		DoubleAngle:     0,                // Depth 3
		K:               16,               // With EphemeralSecretWeight = 32 and 2^{15} slots, ensures < 2^{-138.7} failure probability
		LogMessageRatio: 1,               // q/|m| = 2^10
		Mod1InvDegree:   0,                // Depth 0
	}

	// SlotsToCoeffs parameters (homomorphic decoding)
	SlotsToCoeffsParameters := dft.MatrixLiteral{
		Type:         dft.HomomorphicDecode,
		LogSlots:     params.LogMaxSlots(),
		LogBSGSRatio: 1,
		LevelP:       params.MaxLevelP(),
		Levels:       []int{1, 1}, // qiSlotsToCoeffs
	}

	SlotsToCoeffsParameters.LevelQ = len(SlotsToCoeffsParameters.Levels)

	// Custom bootstrapping.Parameters.
	// All fields are public and can be manually instantiated.
	btpParams := bootstrapping.Parameters{
		ResidualParameters:      params,
		BootstrappingParameters: params,
		SlotsToCoeffsParameters: SlotsToCoeffsParameters,
		Mod1ParametersLiteral:   Mod1ParametersLiteral,
		CoeffsToSlotsParameters: CoeffsToSlotsParameters,
		EphemeralSecretWeight:   32, // > 128bit secure for LogN=16 and LogQP = 115.
		CircuitOrder:            bootstrapping.DecodeThenModUp,
	}


	// We print some information about the residual parameters.
	fmt.Printf("Residual parameters: logN=%d, logSlots=%d, H=%d, sigma=%f, logQP=%f, levels=%d, scale=2^%d\n",
		btpParams.ResidualParameters.LogN(),
		btpParams.ResidualParameters.LogMaxSlots(),
		btpParams.ResidualParameters.XsHammingWeight(),
		btpParams.ResidualParameters.Xe(), params.LogQP(),
		btpParams.ResidualParameters.MaxLevel(),
		btpParams.ResidualParameters.LogDefaultScale())

	// And some information about the bootstrapping parameters.
	// We can notably check that the LogQP of the bootstrapping parameters is smaller than 1550, which ensures
	// 128-bit of security as explained above.
	fmt.Printf("Bootstrapping parameters: logN=%d, logSlots=%d, H(%d; %d), sigma=%f, logQP=%f, levels=%d, scale=2^%d\n",
		btpParams.BootstrappingParameters.LogN(),
		btpParams.BootstrappingParameters.LogMaxSlots(),
		btpParams.BootstrappingParameters.XsHammingWeight(),
		btpParams.EphemeralSecretWeight,
		btpParams.BootstrappingParameters.Xe(),
		btpParams.BootstrappingParameters.LogQP(),
		btpParams.BootstrappingParameters.QCount(),
		btpParams.BootstrappingParameters.LogDefaultScale())

	// Scheme context and keys

	//===========================
	//=== 4) KEYGEN & ENCRYPT ===
	//===========================

	// Now that both the residual and bootstrapping parameters are instantiated, we can
	// instantiate the usual necessary object to encode, encrypt and decrypt.

	// Scheme context and keys
	kgen := rlwe.NewKeyGenerator(params)

	sk, pk := kgen.GenKeyPairNew()

	encoder := ckks.NewEncoder(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	encryptor := rlwe.NewEncryptor(params, pk)

	fmt.Println()
	fmt.Println("Generating bootstrapping evaluation keys...")
	evk, _, err := btpParams.GenEvaluationKeys(sk)
	if err != nil {
		panic(err)
	}
	fmt.Println("Done")

	//========================
	//=== 5) BOOTSTRAPPING ===
	//========================
	iv := make([]uint8, 16)
	symmetricKey := []byte{
        0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00,
    }
	aes, _ := ckks_cipher.NewAESCtr(symmetricKey, params, btpParams, evk, encoder, encryptor, decryptor, iv) 
	// Start1 := time.Now()
	// aes.DebugTest(symmetricKey, 128 )
	aes.HEDecrypt(symmetricKey, 128 )
	// aes.DebugPrint(ctxt[0], "\n\n Finished")
	// end1 := time.Since(Start1)
	// fmt.Printf("The Whole Decryption took %v\n", end1)
}

